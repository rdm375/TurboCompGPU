#!/usr/bin/env python3
"""Controlled leave-one-speed-line-out ablation of compressor-map coordinates."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

TARGETS = ("outlet_total_temperature", "outlet_total_pressure", "specific_work")
MODES = ("mass_coordinate", "boundary_coordinate")


class MLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(2, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 3),
        )

    def forward(self, x):
        return self.layers(x)


def f(row, key):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Bad or missing {key}: {row.get(key)!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {key}")
    return value


def metric(actual, pred):
    rel = (pred - actual) / actual
    return {
        "max_relative_error": float(np.max(np.abs(rel))),
        "p99_relative_error": float(np.percentile(np.abs(rel), 99)),
        "rms_relative_error": float(np.sqrt(np.mean(rel ** 2))),
        "rmse": float(np.sqrt(np.mean((pred - actual) ** 2))),
    }


def fit_predict(x_train, y_train, x_test, *, epochs, lr, seed, device):
    # Match initialization, optimizer, and fold-specific target scaling across modes.
    torch.manual_seed(seed)
    log_y = np.log(y_train).astype(np.float32)
    mean = log_y.mean(axis=0)
    std = log_y.std(axis=0)
    std[std < 1e-8] = 1.0
    y_scaled = (log_y - mean) / std

    model = MLP().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()
    xt = torch.as_tensor(x_train, dtype=torch.float32, device=device)
    yt = torch.as_tensor(y_scaled, dtype=torch.float32, device=device)
    for _ in range(epochs):
        optimizer.zero_grad(set_to_none=True)
        loss = loss_fn(model(xt), yt)
        loss.backward()
        optimizer.step()
    model.eval()
    with torch.no_grad():
        xp = torch.as_tensor(x_test, dtype=torch.float32, device=device)
        log_pred = model(xp).cpu().numpy().astype(np.float64)
    return np.exp(log_pred * std + mean)


def read_dataset(path, epochs, lr, seed, threads):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No rows in {path}")
    required = {
        "model", "rotational_speed_rpm", "mass_flow", "flow_fraction",
        "zero_recovery_flow_kg_s", "critical_diffuser_capacity_kg_s",
        "converged", "validated", "basic_physics_valid", *TARGETS,
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    for row in rows:
        for flag in ("converged", "validated", "basic_physics_valid"):
            if row[flag].strip().lower() not in {"1", "true", "yes"}:
                raise ValueError(f"Invalid sample: {flag} is false")
        if row.get("stage_status", "Complete") != "Complete":
            raise ValueError(f"Invalid stage status: {row.get('stage_status')}")

    rpm_min = min(f(r, "rotational_speed_rpm") for r in rows)
    rpm_max = max(f(r, "rotational_speed_rpm") for r in rows)
    if rpm_max <= rpm_min:
        raise ValueError("Need multiple speed lines")
    for row in rows:
        rpm = f(row, "rotational_speed_rpm")
        row["_speed_x"] = 2 * (rpm - rpm_min) / (rpm_max - rpm_min) - 1
        lo = f(row, "zero_recovery_flow_kg_s")
        hi = f(row, "critical_diffuser_capacity_kg_s")
        mass = f(row, "mass_flow")
        if hi <= lo:
            raise ValueError("Invalid lower/upper flow boundaries")
        frac = (mass - lo) / (hi - lo)
        supplied = f(row, "flow_fraction")
        if abs(frac - supplied) > 2e-8 or frac < -2e-8 or frac > 1 + 2e-8:
            raise ValueError(
                f"Boundary fraction mismatch at {row['model']}, {rpm:g} rpm: "
                f"computed={frac:.12g}, supplied={supplied:.12g}"
            )
        row["_fraction"] = min(1.0, max(0.0, frac))

    torch.set_num_threads(threads)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    report = {
        "experiment": "paired input-coordinate ablation",
        "controlled_settings": {
            "architecture": "2->64->64->32->3 Linear/ReLU hidden",
            "target_transform": "standardized log of positive outputs",
            "optimizer": "Adam",
            "learning_rate": lr,
            "epochs": epochs,
            "cross_validation": "leave one complete speed line out",
            "seed_rule": "same fold seed for both coordinate modes",
            "device": str(device),
        },
        "coordinate_modes": {
            "mass_coordinate": "normalized speed and globally normalized mass flow per gas model",
            "boundary_coordinate": "normalized speed and normalized flow fraction from the two boundary curves",
        },
        "speed_domain_rpm": [rpm_min, rpm_max],
        "models": {},
    }
    predictions = {}

    for model_name in sorted({r["model"] for r in rows}):
        group = [r for r in rows if r["model"] == model_name]
        mass_min = min(f(r, "mass_flow") for r in group)
        mass_max = max(f(r, "mass_flow") for r in group)
        if mass_max <= mass_min:
            raise ValueError(f"Need multiple mass-flow points for {model_name}")
        x_modes = {
            "mass_coordinate": np.array([
                [r["_speed_x"], 2 * (f(r, "mass_flow") - mass_min) / (mass_max - mass_min) - 1]
                for r in group
            ], dtype=np.float32),
            "boundary_coordinate": np.array([
                [r["_speed_x"], 2 * r["_fraction"] - 1]
                for r in group
            ], dtype=np.float32),
        }
        y = np.array([[f(r, t) for t in TARGETS] for r in group], dtype=np.float64)
        speed_ids = np.array([int(round(f(r, "rotational_speed_rpm"))) for r in group])
        speeds = sorted(set(speed_ids.tolist()))
        if len(speeds) < 3:
            raise ValueError(f"Need at least 3 speed lines for {model_name}")

        report["models"][model_name] = {"samples": len(group), "speed_lines": len(speeds), "modes": {}}
        for mode_index, mode in enumerate(MODES):
            pred = np.empty_like(y)
            for fold, held_speed in enumerate(speeds):
                train_mask = speed_ids != held_speed
                test_mask = ~train_mask
                pred[test_mask] = fit_predict(
                    x_modes[mode][train_mask], y[train_mask], x_modes[mode][test_mask],
                    epochs=epochs, lr=lr, seed=seed + fold,
                    device=device,
                )

            mode_metrics = {t: metric(y[:, j], pred[:, j]) for j, t in enumerate(TARGETS)}
            endpoints = {}
            fractions = np.array([r["_fraction"] for r in group])
            for boundary, label in ((0.0, "zero_recovery"), (1.0, "critical_capacity")):
                mask = np.isclose(fractions, boundary, atol=1e-7, rtol=0.0)
                endpoints[label] = {
                    t: metric(y[mask, j], pred[mask, j])
                    for j, t in enumerate(TARGETS)
                }
            report["models"][model_name]["modes"][mode] = {
                "mass_flow_normalization": [mass_min, mass_max],
                "metrics": mode_metrics,
                "boundary_endpoint_metrics": endpoints,
            }
            predictions[(model_name, mode)] = (group, pred)

    args_out = Path(args.output_dir)
    args_out.mkdir(parents=True, exist_ok=True)
    for (model_name, mode), (group, pred) in predictions.items():
        outdir = args_out / mode
        outdir.mkdir(parents=True, exist_ok=True)
        with (outdir / f"{model_name.lower()}_holdout.csv").open("w", newline="") as stream:
            fields = ["model", "held_out_speed_rpm", "mass_flow", "flow_fraction"]
            for t in TARGETS:
                fields.extend((f"{t}_actual", f"{t}_predicted", f"{t}_relative_error"))
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row, predicted in zip(group, pred):
                output = {
                    "model": model_name,
                    "held_out_speed_rpm": f(row, "rotational_speed_rpm"),
                    "mass_flow": f(row, "mass_flow"),
                    "flow_fraction": row["_fraction"],
                }
                for j, t in enumerate(TARGETS):
                    actual = f(row, t)
                    output[f"{t}_actual"] = actual
                    output[f"{t}_predicted"] = predicted[j]
                    output[f"{t}_relative_error"] = (predicted[j] - actual) / actual
                writer.writerow(output)

    # Compact per-target deltas: negative means the boundary coordinates reduced error.
    report["max_error_delta_boundary_minus_mass"] = {}
    for model_name, model_data in report["models"].items():
        report["max_error_delta_boundary_minus_mass"][model_name] = {}
        for target in TARGETS:
            m = model_data["modes"]["mass_coordinate"]["metrics"][target]["max_relative_error"]
            b = model_data["modes"]["boundary_coordinate"]["metrics"][target]["max_relative_error"]
            report["max_error_delta_boundary_minus_mass"][model_name][target] = b - m

    report_path = args_out / "coordinate_ablation_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Device: {device}; folds: leave one complete speed line out")
    for model_name, model_data in report["models"].items():
        print(f"{model_name}")
        for target in TARGETS:
            vals = []
            for mode in MODES:
                m = model_data["modes"][mode]["metrics"][target]
                vals.append(f"{mode} max={100*m['max_relative_error']:.4f}% "
                            f"rms={100*m['rms_relative_error']:.4f}%")
            print(f"  {target}: " + " | ".join(vals))
    print(f"Report: {report_path}")
    print(f"Holdout CSVs: {args_out}/{{mass_coordinate,boundary_coordinate}}/")
    return 0


if __name__ == "__main__":
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
    raise SystemExit(read_dataset(args.csv, args.epochs, args.learning_rate,
                                 args.seed, args.cpu_threads))
