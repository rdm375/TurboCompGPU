#!/usr/bin/env python3
"""Fit a ReLU surrogate to TurboCompGPU's existing compressor-map samples.

This is a first, fixed-condition interpolation experiment. It trains a
separate network per thermodynamic model from mass flow and speed to total
outlet temperature, total outlet pressure, and specific work.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from torch import nn


FEATURES = ("mass_flow", "rotational_speed_rpm")
TARGETS = (
    "outlet_total_temperature",
    "outlet_total_pressure",
    "specific_work",
)


class FlowEnvelope:
    """Piecewise-linear lower/upper mass-flow limits as functions of speed."""

    def __init__(self, speed: np.ndarray, lower: np.ndarray,
                 upper: np.ndarray, source: str) -> None:
        order = np.argsort(speed)
        self.speed = np.asarray(speed, dtype=np.float64)[order]
        self.lower = np.asarray(lower, dtype=np.float64)[order]
        self.upper = np.asarray(upper, dtype=np.float64)[order]
        self.source = source
        if len(self.speed) < 2 or np.any(np.diff(self.speed) <= 0.0):
            raise ValueError("flow envelope needs at least two unique speeds")
        if (not np.all(np.isfinite(self.speed)) or
                not np.all(np.isfinite(self.lower)) or
                not np.all(np.isfinite(self.upper)) or
                np.any(self.lower <= 0.0) or
                np.any(self.upper <= self.lower)):
            raise ValueError("flow envelope contains invalid limits")

    def limits(self, speed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        speed = np.asarray(speed, dtype=np.float64)
        tolerance = 1.0e-10 * max(1.0, abs(self.speed[-1]))
        if (np.any(speed < self.speed[0] - tolerance) or
                np.any(speed > self.speed[-1] + tolerance)):
            raise ValueError(
                "requested speed is outside the supplied flow-envelope range")
        return (np.interp(speed, self.speed, self.lower),
                np.interp(speed, self.speed, self.upper))

    def flow_fraction(self, x: np.ndarray) -> np.ndarray:
        lower, upper = self.limits(x[:, 1])
        fraction = (x[:, 0] - lower) / (upper - lower)
        tolerance = 1.0e-8
        if np.any(fraction < -tolerance) or np.any(fraction > 1.0 + tolerance):
            raise ValueError("sample lies outside the configured flow envelope")
        return np.clip(fraction, 0.0, 1.0)


def read_flow_envelope(path: Path) -> FlowEnvelope:
    required = {"rotational_speed_rpm", "mass_flow_min", "mass_flow_max"}
    rows = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            try:
                rows.append((float(row["rotational_speed_rpm"]),
                             float(row["mass_flow_min"]),
                             float(row["mass_flow_max"])))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_no}: invalid envelope row") from exc
    if not rows:
        raise ValueError(f"{path}: no envelope rows")
    values = np.asarray(rows, dtype=np.float64)
    if len(np.unique(values[:, 0])) != len(values):
        raise ValueError(f"{path}: speed values must be unique")
    return FlowEnvelope(values[:, 0], values[:, 1], values[:, 2],
                        f"user-supplied:{path}")


def write_flow_envelope(path: Path, envelope: FlowEnvelope) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("rotational_speed_rpm", "mass_flow_min",
                         "mass_flow_max"))
        writer.writerows(zip(envelope.speed, envelope.lower, envelope.upper))


def encoded_features(x: np.ndarray, envelope: FlowEnvelope,
                     speed_domain: tuple[float, float]) -> np.ndarray:
    """Map physical coordinates to flow fraction and normalized speed."""
    speed_lo, speed_hi = speed_domain
    tolerance = 1.0e-10 * max(1.0, abs(speed_hi))
    if (np.any(x[:, 1] < speed_lo - tolerance) or
            np.any(x[:, 1] > speed_hi + tolerance)):
        raise ValueError("requested speed is outside the trained domain")
    flow_fraction = envelope.flow_fraction(x)
    speed_fraction = (x[:, 1] - speed_lo) / (speed_hi - speed_lo)
    return np.column_stack((flow_fraction, speed_fraction)).astype(np.float32)


class ReluMap(nn.Module):
    def __init__(self, width: int) -> None:
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(2, width),
            nn.ReLU(),
            nn.Linear(width, width),
            nn.ReLU(),
            nn.Linear(width, len(TARGETS)),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.layers(x)


def read_samples(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: dict[str, dict[tuple[float, float], np.ndarray]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"model", *FEATURES, *TARGETS,
                    "converged", "validated", "basic_physics_valid"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")

        for line_no, row in enumerate(reader, start=2):
            valid = all(row[name].strip().lower() in {"1", "true"}
                        for name in ("converged", "validated",
                                     "basic_physics_valid"))
            if not valid:
                raise ValueError(f"{path}:{line_no}: sample did not pass "
                                 "solver validation")

            model = row["model"].strip()
            x = tuple(float(row[name]) for name in FEATURES)
            y = np.array([float(row[name]) for name in TARGETS],
                         dtype=np.float64)
            if not np.all(np.isfinite((*x, *y))) or np.any(y <= 0.0):
                raise ValueError(f"{path}:{line_no}: non-finite or "
                                 "non-positive feature/target")

            key = (x[0], x[1])
            previous = grouped[model].get(key)
            if previous is not None:
                if not np.allclose(previous, y, rtol=1.0e-8, atol=1.0e-8):
                    raise ValueError(
                        f"{path}:{line_no}: repeated condition {key} "
                        "has inconsistent targets")
                # Independent and continuation rows represent the same
                # operating condition. Keep one sample, not two.
                continue
            grouped[model][key] = y

    result = {}
    for model, points in grouped.items():
        ordered = sorted(points.items())
        x = np.asarray([point[0] for point in ordered], dtype=np.float64)
        y = np.asarray([point[1] for point in ordered], dtype=np.float64)
        if len(np.unique(x[:, 1])) < 3:
            raise ValueError(f"{model}: need at least three speed lines")
        result[model] = (x, y)
    if not result:
        raise ValueError(f"No valid samples found in {path}")
    return result


def fit_model(
    x: np.ndarray,
    y: np.ndarray,
    envelope: FlowEnvelope,
    speed_domain: tuple[float, float],
    width: int,
    epochs: int,
    learning_rate: float,
    seed: int,
    cpu_threads: int,
) -> tuple[ReluMap, np.ndarray, np.ndarray]:
    torch.manual_seed(seed)
    torch.set_num_threads(cpu_threads)

    # The first coordinate follows the compressor envelope, while the second
    # uses the declared speed domain. These bounds remain fixed in every CV
    # fold, including when an edge speed line is held out.
    x_scaled = 2.0 * encoded_features(x, envelope, speed_domain) - 1.0

    # Positive compressor outputs are modeled in log space so relative
    # output errors receive comparable weight and predictions stay positive.
    log_y = np.log(y)
    y_mean = log_y.mean(axis=0)
    y_std = log_y.std(axis=0)
    y_std = np.where(y_std > 1.0e-12, y_std, 1.0)
    y_scaled = ((log_y - y_mean) / y_std).astype(np.float32)

    model = ReluMap(width)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_fn = nn.MSELoss()
    tx = torch.from_numpy(x_scaled)
    ty = torch.from_numpy(y_scaled)

    model.train()
    for _ in range(epochs):
        optimizer.zero_grad(set_to_none=True)
        loss = loss_fn(model(tx), ty)
        loss.backward()
        optimizer.step()

    model.eval()
    return model, y_mean, y_std


def predict(
    model: ReluMap,
    x: np.ndarray,
    envelope: FlowEnvelope,
    speed_domain: tuple[float, float],
    y_mean: np.ndarray,
    y_std: np.ndarray,
) -> np.ndarray:
    with torch.no_grad():
        features = 2.0 * encoded_features(x, envelope, speed_domain) - 1.0
        scaled = model(torch.from_numpy(features))
        log_y = scaled.numpy() * y_std + y_mean
    return np.exp(log_y)


def error_summary(actual: np.ndarray,
                  predicted: np.ndarray) -> dict[str, float]:
    relative = (predicted - actual) / actual
    absolute = predicted - actual
    return {
        "max_relative_error": float(np.max(np.abs(relative))),
        "p99_relative_error": float(np.quantile(np.abs(relative), 0.99)),
        "rms_relative_error": float(np.sqrt(np.mean(relative ** 2))),
        "rmse": float(np.sqrt(np.mean(absolute ** 2))),
    }


def evaluate_speedline_folds(
    x: np.ndarray,
    y: np.ndarray,
    envelope: FlowEnvelope,
    speed_domain: tuple[float, float],
    width: int,
    epochs: int,
    learning_rate: float,
    seed: int,
    cpu_threads: int,
) -> list[dict[str, object]]:
    fold_records: list[dict[str, object]] = []

    speeds = np.unique(x[:, 1])
    for fold_index, speed in enumerate(speeds):
        held = np.isclose(x[:, 1], speed)
        train = ~held
        model, y_mean, y_std = fit_model(
            x[train], y[train], envelope, speed_domain,
            width, epochs, learning_rate,
            seed + fold_index, cpu_threads)
        predicted = predict(model, x[held], envelope, speed_domain,
                            y_mean, y_std)

        for input_row, actual, estimate in zip(
                x[held], y[held], predicted):
            record: dict[str, object] = {
                "held_out_speed_rpm": float(speed),
                "mass_flow": float(input_row[0]),
                "rotational_speed_rpm": float(input_row[1]),
                "flow_fraction": float(envelope.flow_fraction(
                    input_row.reshape(1, -1))[0]),
            }
            for j, target in enumerate(TARGETS):
                record[f"{target}_actual"] = float(actual[j])
                record[f"{target}_predicted"] = float(estimate[j])
                record[f"{target}_relative_error"] = float(
                    (estimate[j] - actual[j]) / actual[j])
            fold_records.append(record)

    return fold_records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        nargs="?",
        type=Path,
        default=Path("benchmarks/results-compressor-map/map.csv"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/relu-surrogate"))
    parser.add_argument(
        "--envelope-csv", type=Path,
        help=("CSV with rotational_speed_rpm,mass_flow_min,mass_flow_max. "
              "Limits are linearly interpolated with speed."))
    parser.add_argument("--width", type=int, default=16)
    parser.add_argument("--epochs", type=int, default=2000)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=375)
    parser.add_argument("--cpu-threads", type=int, default=1)
    args = parser.parse_args()

    if args.width < 2 or args.epochs < 1 or args.cpu_threads < 1:
        parser.error("width, epochs, and cpu-threads must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    by_model = read_samples(args.csv)
    all_x = np.concatenate([values[0] for values in by_model.values()])
    speed_domain = (float(all_x[:, 1].min()), float(all_x[:, 1].max()))
    if args.envelope_csv:
        envelope = read_flow_envelope(args.envelope_csv)
    else:
        # The present benchmark samples a rectangle. Use its known flow range
        # as a provisional support envelope; this is not a surge/choke model.
        lo, hi = float(all_x[:, 0].min()), float(all_x[:, 0].max())
        envelope = FlowEnvelope(
            np.asarray(speed_domain), np.asarray([lo, lo]),
            np.asarray([hi, hi]), "sample-support rectangle (not physical limits)")
    envelope.limits(np.asarray(speed_domain))
    used_envelope_path = args.output_dir / "flow_envelope_used.csv"
    write_flow_envelope(used_envelope_path, envelope)

    all_summary: dict[str, object] = {}
    cv_path = args.output_dir / "speedline_holdout.csv"

    with cv_path.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "model", "held_out_speed_rpm", "mass_flow",
            "rotational_speed_rpm", "flow_fraction",
        ]
        for target in TARGETS:
            fields.extend((f"{target}_actual", f"{target}_predicted",
                           f"{target}_relative_error"))
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()

        for model_index, (model_name, (x, y)) in enumerate(
                sorted(by_model.items())):
            records = evaluate_speedline_folds(
                x, y, envelope, speed_domain,
                args.width, args.epochs, args.learning_rate,
                args.seed + 100 * model_index, args.cpu_threads)
            for record in records:
                writer.writerow({"model": model_name, **record})

            fold_metrics = {}
            for target_index, target in enumerate(TARGETS):
                actual = np.array(
                    [float(r[f"{target}_actual"]) for r in records])
                predicted = np.array(
                    [float(r[f"{target}_predicted"]) for r in records])
                fold_metrics[target] = error_summary(actual, predicted)

            # Fit deployment candidate on all available conditions after
            # completing the structured held-out evaluation.
            final_model, y_mean, y_std = fit_model(
                x, y, envelope, speed_domain,
                args.width, args.epochs, args.learning_rate,
                args.seed + 1000 + model_index, args.cpu_threads)
            checkpoint = {
                "model_name": model_name,
                "architecture": [2, args.width, args.width, len(TARGETS)],
                "features": list(FEATURES),
                "targets": list(TARGETS),
                "state_dict": final_model.state_dict(),
                "feature_encoding": [
                    "flow_fraction=(mass_flow-min_flow(speed))/(max_flow(speed)-min_flow(speed))",
                    "speed_fraction=(speed-speed_min)/(speed_max-speed_min)",
                ],
                "speed_domain_rpm": np.asarray(speed_domain),
                "envelope_speed_rpm": envelope.speed,
                "envelope_mass_flow_min": envelope.lower,
                "envelope_mass_flow_max": envelope.upper,
                "envelope_source": envelope.source,
                "log_target_mean": y_mean,
                "log_target_std": y_std,
            }
            safe_model_name = {
                "IdealGas": "ideal_gas",
                "PengRobinson": "peng_robinson",
            }.get(model_name, model_name.lower().replace(" ", "_"))
            torch.save(
                checkpoint,
                args.output_dir / f"{safe_model_name}.pt")

            all_summary[model_name] = {
                "unique_operating_points": int(len(x)),
                "feature_ranges": {
                    feature: [float(x[:, i].min()), float(x[:, i].max())]
                    for i, feature in enumerate(FEATURES)
                },
                "flow_envelope": {
                    "source": envelope.source,
                    "profile_csv": str(used_envelope_path),
                    "interpolation": "piecewise linear in rotational speed",
                    "limits_by_speed": [
                        {"rotational_speed_rpm": float(n),
                         "mass_flow_min": float(mn),
                         "mass_flow_max": float(mx)}
                        for n, mn, mx in zip(
                            envelope.speed, envelope.lower, envelope.upper)
                    ],
                },
                "structured_validation": {
                    "method": "leave one rotational-speed line out",
                    "folds": int(len(np.unique(x[:, 1]))),
                    "note": (
                        "The lowest and highest speed folds are extrapolation "
                        "checks; interior speed folds test interpolation."
                    ),
                    "metrics_by_target": fold_metrics,
                },
                "checkpoint": f"{safe_model_name}.pt",
            }

    report = {
        "label_source": str(args.csv),
        "surrogate_type": "fully connected ReLU network",
        "architecture": f"2->{args.width}->{args.width}->3",
        "features": list(FEATURES),
        "targets": list(TARGETS),
        "target_transform": "natural log; exponentiated at inference",
        "dataset_note": (
            "Rows with identical model, mass flow, and speed are deduplicated "
            "after verifying that their target values agree. For this source, "
            "inlet conditions, geometry, slip, prewhirl, and entropy generation "
            "are fixed by the benchmark and are not surrogate inputs."
        ),
        "input_encoding": {
            "features": ["flow fraction within envelope", "speed fraction"],
            "out_of_envelope": "rejected during feature encoding",
            "speed_domain_rpm": list(speed_domain),
        },
        "validation": "leave-one-speed-line-out cross-validation",
        "models": all_summary,
        "qualification": (
            "Preliminary interpolation only. Samples are generated by the "
            "current numerical model, not the legacy hackathon simulator or "
            "experimental compressor measurements."
        ),
    }
    report_path = args.output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")

    print(f"Samples: {args.csv}")
    for model_name, summary in all_summary.items():
        print(f"{model_name}: "
              f"{summary['unique_operating_points']} unique conditions")
        for target, metrics in summary["structured_validation"][
                "metrics_by_target"].items():
            print(
                f"  {target}: max_rel="
                f"{metrics['max_relative_error']:.6g} "
                f"p99_rel={metrics['p99_relative_error']:.6g} "
                f"rms_rel={metrics['rms_relative_error']:.6g}"
            )
        print(f"  checkpoint={summary['checkpoint']}")
    print(f"Cross-validation predictions: {cv_path}")
    print(f"Flow envelope profile: {used_envelope_path} ({envelope.source})")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
