#!/usr/bin/env python3
"""Fit an effective head parameter on Eckardt B training lines and gate a holdout.

This is a structured speed-line holdout, not a blind compressor-family test.
The fitted value is an empirical effective slip/head parameter and must not be
interpreted as a measured physical slip factor.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path


def read_points(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    return list(csv.DictReader(lines, fieldnames=[
        "rpm", "mass_flow_kg_s", "measured_total_pressure_ratio"
    ]))


def write_points(path: Path, points: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        for point in points:
            writer.writerow((point["rpm"], point["mass_flow_kg_s"],
                             point["measured_total_pressure_ratio"]))


def evaluate(executable: Path, points: list[dict[str, str]], slip: float,
             scratch: Path) -> list[dict[str, str]]:
    source = scratch / "input.csv"
    output = scratch / "prediction.csv"
    write_points(source, points)
    run = subprocess.run(
        [str(executable), str(source), str(output), f"{slip:.12g}", "B"],
        check=False, text=True, capture_output=True,
    )
    if run.returncode:
        raise RuntimeError(
            f"Eckardt B evaluator failed (exit {run.returncode}):\n"
            f"{run.stdout}{run.stderr}"
        )
    with output.open(newline="") as stream:
        return list(csv.DictReader(stream))


def errors(rows: list[dict[str, str]]) -> list[float]:
    return [float(row["relative_error_percent"]) for row in rows
            if row["stage_complete"].lower() in {"1", "true"}]


def objective(executable: Path, train: list[dict[str, str]], slip: float,
              scratch: Path) -> float:
    values = errors(evaluate(executable, train, slip, scratch))
    if len(values) != len(train):
        return math.inf
    return sum(value * value for value in values) / len(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--map", type=Path,
                        default=Path("data/validation/eckardt_impeller_b.csv"))
    parser.add_argument("--train-rpm", nargs="+", type=int,
                        default=[10000, 12000, 14000])
    parser.add_argument("--holdout-rpm", nargs="+", type=int, default=[16000])
    parser.add_argument("--max-holdout-error-percent", type=float, default=7.0)
    parser.add_argument("--slip-bounds", nargs=2, type=float,
                        default=[0.60, 0.90], metavar=("LOW", "HIGH"))
    parser.add_argument("--output", type=Path,
                        default=Path("benchmarks/results-eckardt-impeller-b/"
                                     "speedline_holdout.csv"))
    parser.add_argument("--report", type=Path,
                        default=Path("benchmarks/results-eckardt-impeller-b/"
                                     "speedline_holdout.md"))
    args = parser.parse_args()

    all_points = read_points(args.map)
    training = [p for p in all_points if int(float(p["rpm"])) in args.train_rpm]
    heldout = [p for p in all_points if int(float(p["rpm"])) in args.holdout_rpm]
    if not training or not heldout:
        raise SystemExit("Training or holdout speed lines are empty")
    if set(args.train_rpm) & set(args.holdout_rpm):
        raise SystemExit("Training and holdout speed lines must be disjoint")
    low, high = args.slip_bounds
    if not (0.0 < low < high <= 1.0):
        raise SystemExit("Slip bounds must satisfy 0 < LOW < HIGH <= 1")

    # Deterministic bounded golden-section minimization of training-line mean
    # squared relative pressure-ratio error. Holdout rows are not evaluated
    # until the fitted parameter has been selected.
    with tempfile.TemporaryDirectory(prefix="eckardt-b-holdout-") as temp:
        scratch = Path(temp)
        ratio = (math.sqrt(5.0) - 1.0) / 2.0
        left, right = low, high
        x1, x2 = right - ratio * (right - left), left + ratio * (right - left)
        f1 = objective(args.executable, training, x1, scratch)
        f2 = objective(args.executable, training, x2, scratch)
        for _ in range(32):
            if f1 <= f2:
                right, x2, f2 = x2, x1, f1
                x1 = right - ratio * (right - left)
                f1 = objective(args.executable, training, x1, scratch)
            else:
                left, x1, f1 = x1, x2, f2
                x2 = left + ratio * (right - left)
                f2 = objective(args.executable, training, x2, scratch)
        fitted = 0.5 * (left + right)
        train_rows = evaluate(args.executable, training, fitted, scratch)
        holdout_rows = evaluate(args.executable, heldout, fitted, scratch)
        all_rows = evaluate(args.executable, all_points, fitted, scratch)

    holdout_errors = errors(holdout_rows)
    train_errors = errors(train_rows)
    if len(holdout_errors) != len(heldout) or len(train_errors) != len(training):
        raise SystemExit("One or more training/holdout points did not converge")
    holdout_max = max(abs(value) for value in holdout_errors)
    holdout_mape = sum(abs(value) for value in holdout_errors) / len(holdout_errors)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]) + ["split"])
        writer.writeheader()
        for row in all_rows:
            row["split"] = (
                "training" if int(float(row["rpm"])) in args.train_rpm
                else "held_out" if int(float(row["rpm"])) in args.holdout_rpm
                else "not_scored"
            )
            writer.writerow(row)

    passed = holdout_max <= args.max_holdout_error_percent
    report = [
        "# Eckardt B speed-line holdout",
        "",
        "**Scope:** this is a same-family, structured speed-line holdout. The",
        "effective slip/head parameter is fitted on the listed training lines.",
        "It is empirical and must not be interpreted as a physical slip",
        "measurement or a cross-family qualification.",
        "",
        f"- Training speed lines: {', '.join(map(str, args.train_rpm))} rpm ({len(training)} points)",
        f"- Held-out speed lines: {', '.join(map(str, args.holdout_rpm))} rpm ({len(heldout)} points)",
        f"- Fitted effective parameter: {fitted:.8f}",
        f"- Training MAPE: {sum(abs(x) for x in train_errors) / len(train_errors):.3f}%",
        f"- Held-out MAPE: {holdout_mape:.3f}%",
        f"- Held-out maximum absolute error: {holdout_max:.3f}%",
        f"- Predeclared held-out maximum-error gate: {args.max_holdout_error_percent:.3f}%",
        f"- Result: **{'PASS' if passed else 'FAIL'}**",
        "",
        "The 7% default is a provisional project gate informed by the maximum",
        "Eckardt B pressure-ratio deviation reported for a separate published",
        "mean-line model. It is not an experimental uncertainty bound. Change",
        "the gate only before fitting, and record the reason.",
        "",
        f"Pointwise results: `{args.output}`.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
