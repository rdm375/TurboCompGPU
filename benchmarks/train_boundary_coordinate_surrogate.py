#!/usr/bin/env python3
"""Train a complete-stage surrogate in normalized speed/operating-envelope coordinates."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

TARGETS = (
    "outlet_total_temperature",
    "outlet_total_pressure",
    "specific_work",
)


class BoundaryCoordinateMLP(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(2, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 3),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.layers(x)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    required = {
        "model", "rotational_speed_rpm", "mass_flow", "flow_fraction",
        "zero_recovery_flow_kg_s", "critical_diffuser_capacity_kg_s",
        "converged", "validated", "basic_physics_valid",
        *TARGETS,
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return rows


def as_float(row: dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (TypeError, ValueError, KeyError) as exc:
        raise ValueError(f"Invalid {key} value {row.get(key)!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {key} value")
    return value


def require_true(row: dict[str, str], key: str) -> bool:
    return row[key].strip().lower() in {"1", "true", "yes"}


def prepare(rows: list[dict[str, str]], speed_min: float, speed_max: float):
    """Return arrays with x=[normalized speed, normalized envelope fraction]."""
    x_values, y_values, metadata = [], [], []
    group_bounds: dict[tuple[str, int], tuple[float, float]] = {}
    seen = set()

    for row in rows:
        for flag in ("converged", "validated", "basic_physics_valid"):
            if not require_true(row, flag):
                raise ValueError(f"Rejected {row['model']} row: {flag}=false")
        if row.get("stage_status", "Complete") != "Complete":
            raise ValueError(f"Rejected row with stage_status={row['stage_status']}")

        model = row["model"]
        rpm = as_float(row, "rotational_speed_rpm")
        flow = as_float(row, "mass_flow")
        lo = as_float(row, "zero_recovery_flow_kg_s")
        hi = as_float(row, "critical_diffuser_capacity_kg_s")
        if not hi > lo:
            raise ValueError(f"Invalid flow bounds at {model}, {rpm:g} rpm")

        speed_id = int(round(rpm))
        key = (model, speed_id)
        old = group_bounds.setdefault(key, (lo, hi))
        if not (math.isclose(old[0], lo, rel_tol=0, abs_tol=1e-9) and
                math.isclose(old[1], hi, rel_tol=0, abs_tol=1e-9)):
            raise ValueError(f"Inconsistent boundary values for {key}")

        fraction = (flow - lo) / (hi - lo)
        supplied_fraction = as_float(row, "flow_fraction")
        if not math.isclose(fraction, supplied_fraction, rel_tol=0, abs_tol=2e-8):
            raise ValueError(
                f"flow_fraction disagrees with mass-flow bounds at {key}: "
                f"computed={fraction:.12g}, supplied={supplied_fraction:.12g}"
            )
        if fraction < -2e-8 or fraction > 1.0 + 2e-8:
            raise ValueError(f"Sample lies outside operating envelope: {key}, {fraction}")
        fraction = min(1.0, max(0.0, fraction))

        row_key = (model, speed_id, round(flow, 9))
        if row_key in seen:
            raise ValueError(f"Duplicate operating point: {row_key}")
        seen.add(row_key)

        speed_coordinate = 2.0 * (rpm - speed_min) / (speed_max - speed_min) - 1.0
        x_values.append([speed_coordinate, 2.0 * fraction - 1.0])
        y_values.append([as_float(row, target) for target in TARGETS])
        metadata.append({
            "model": model,
            "rpm": rpm,
            "mass_flow": flow,
            "flow_fraction": fraction,
            "lower": lo,
            "upper": hi,
        })

    return np.asarray(x_values, dtype=np.float32), np.asarray(y_values, dtype=np.float64), metadata


def fit_model(x: np.ndarray, y: np.ndarray, *, epochs: int, seed: int,
              learning_rate: float, device: torch.device, threads: int):
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(threads)

    # Log targets make squared error correspond to relative output error.
    log_y = np.log(y).astype(np.float32)
    target_mean = log_y.mean(axis=0)
    target_std = log_y.std(axis=0)
    target_std[target_std < 1e-8] = 1.0
    normalized_y = (log_y - target_mean) / target_std

    x_tensor = torch.as_tensor(x, dtype=torch.float32, device=device)
    y_tensor = torch.as_tensor(normalized_y, dtype=torch.float32, device=device)
    model = BoundaryCoordinateMLP().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_fn = nn.MSELoss()
    for _ in range(epochs):
        optimizer.zero_grad(set_to_none=True)
        loss = loss_fn(model(x_tensor), y_tensor)
        loss.backward()
        optimizer.step()

    return model, target_mean.astype(np.float64), target_std.astype(np.float64)


def predict(model, x: np.ndarray, target_mean: np.ndarray,
            target_std: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        xt = torch.as_tensor(x, dtype=torch.float32, device=device)
        log_pred = model(xt).cpu().numpy().astype(np.float64)
    return np.exp(log_pred * target_std + target_mean)


def metric_summary(actual: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    rel = (predicted - actual) / actual
    abs_rel = np.abs(rel)
    return {
        "max_relative_error": float(abs_rel.max()),
        "p99_relative_error": float(np.percentile(abs_rel, 99)),
        "rms_relative_error": float(np.sqrt(np.mean(rel * rel))),
        "rmse": float(np.sqrt(np.mean((predicted - actual) ** 2))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--cpu-threads", type=int, default=1)
    args = parser.parse_args()

    if args.epochs < 1 or args.cpu_threads < 1 or args.learning_rate <= 0:
        parser.error("epochs, cpu-threads, and learning-rate must be positive")

    rows = read_rows(args.csv)
    speeds = np.array([as_float(r, "rotational_speed_rpm") for r in rows])
    speed_min, speed_max = float(speeds.min()), float(speeds.max())
    if speed_max <= speed_min:
        raise ValueError("Need at least two distinct speeds")
    x_all, y_all, meta_all = prepare(rows, speed_min, speed_max)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prediction_rows = []
    report = {
        "label_source": "complete-stage solver; supplied boundary-grid CSV",
        "architecture": "2->64->64->32->3 Linear/ReLU hidden; log targets",
        "input_features": [
            "normalized rotational speed in [-1,1]",
            "normalized flow fraction in [-1,1]",
        ],
        "flow_fraction_definition":
            "(mass_flow - zero_recovery_flow)/(critical_diffuser_capacity - zero_recovery_flow)",
        "speed_domain_rpm": [speed_min, speed_max],
        "flow_fraction_domain": [0.0, 1.0],
        "training_device": str(device),
        "epochs": args.epochs,
        "learning_rate": args.learning_rate,
        "cross_validation": "leave one complete speed line out; all rows at each RPM held out together",
        "models": {},
    }

    for model_name in sorted({m["model"] for m in meta_all}):
        idx = np.array([i for i, m in enumerate(meta_all) if m["model"] == model_name])
        model_speeds = sorted({int(round(meta_all[i]["rpm"])) for i in idx})
        if len(model_speeds) < 3:
            raise ValueError(f"Need at least 3 speed lines for {model_name}")

        actual_all = y_all[idx]
        pred_all = np.empty_like(actual_all)
        for fold, held_speed in enumerate(model_speeds):
            held_mask = np.array([int(round(meta_all[i]["rpm"])) == held_speed for i in idx])
            train_idx = idx[~held_mask]
            test_idx = idx[held_mask]
            net, y_mean, y_std = fit_model(
                x_all[train_idx], y_all[train_idx], epochs=args.epochs,
                seed=args.seed + fold, learning_rate=args.learning_rate,
                device=device, threads=args.cpu_threads)
            pred_all[held_mask] = predict(
                net, x_all[test_idx], y_mean, y_std, device)

        model_report = {"samples": int(len(idx)), "speed_lines": len(model_speeds), "metrics": {}}
        for j, target in enumerate(TARGETS):
            model_report["metrics"][target] = metric_summary(actual_all[:, j], pred_all[:, j])

        endpoint_metrics = {}
        model_meta = [meta_all[i] for i in idx]
        fractions = np.array([m["flow_fraction"] for m in model_meta])
        for endpoint, label in ((0.0, "zero_recovery"), (1.0, "critical_capacity")):
            mask = np.isclose(fractions, endpoint, atol=1e-7, rtol=0.0)
            endpoint_metrics[label] = {
                target: metric_summary(actual_all[mask, j], pred_all[mask, j])
                for j, target in enumerate(TARGETS)
            }
        model_report["boundary_endpoint_metrics"] = endpoint_metrics
        report["models"][model_name] = model_report

        for local_i, global_i in enumerate(idx):
            meta = meta_all[global_i]
            row = {
                "model": model_name,
                "held_out_speed_rpm": meta["rpm"],
                "mass_flow": meta["mass_flow"],
                "flow_fraction": meta["flow_fraction"],
            }
            for j, target in enumerate(TARGETS):
                actual = float(actual_all[local_i, j])
                predicted = float(pred_all[local_i, j])
                row[f"{target}_actual"] = actual
                row[f"{target}_predicted"] = predicted
                row[f"{target}_relative_error"] = (predicted - actual) / actual
            prediction_rows.append(row)

        final_model, final_mean, final_std = fit_model(
            x_all[idx], y_all[idx], epochs=args.epochs, seed=args.seed,
            learning_rate=args.learning_rate, device=device, threads=args.cpu_threads)
        torch.save({
            "state_dict": final_model.cpu().state_dict(),
            "target_log_mean": final_mean,
            "target_log_std": final_std,
            "input_features": report["input_features"],
            "speed_domain_rpm": [speed_min, speed_max],
            "flow_fraction_domain": [0.0, 1.0],
            "targets": list(TARGETS),
            "model": model_name,
        }, args.output_dir / f"{model_name.lower()}.pt")

        print(f"{model_name}: {len(idx)} samples, {len(model_speeds)} held-out speed lines")
        for target in TARGETS:
            m = model_report["metrics"][target]
            print(f"  {target}: max_rel={m['max_relative_error']:.6g} "
                  f"p99_rel={m['p99_relative_error']:.6g} "
                  f"rms_rel={m['rms_relative_error']:.6g}")
        print(f"  checkpoint={model_name.lower()}.pt")

    pred_path = args.output_dir / "speedline_holdout.csv"
    with pred_path.open("w", newline="") as stream:
        fieldnames = ["model", "held_out_speed_rpm", "mass_flow", "flow_fraction"]
        for target in TARGETS:
            fieldnames.extend((f"{target}_actual", f"{target}_predicted", f"{target}_relative_error"))
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(prediction_rows)

    report_path = args.output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Cross-validation predictions: {pred_path}")
    print(f"Report: {report_path}")
    print("Boundary coordinate: explicit speed-dependent zero-recovery/capacity bounds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
