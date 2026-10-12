#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

root = (
    Path(__file__).resolve().parent
    / "results-characteristic-convergence"
)


def load(filename):
    with (root / filename).open(newline="") as stream:
        return list(csv.DictReader(stream))


optimizer = load("optimizer_convergence.csv")
derivatives = load("derivative_convergence.csv")

opt_groups = defaultdict(list)
der_groups = defaultdict(list)

for row in optimizer:
    opt_groups[float(row["rpm"])].append(row)

for row in derivatives:
    der_groups[float(row["rpm"])].append(row)

print("\nSTAGE 7H OPTIMIZER CONVERGENCE")
print("=" * 100)

print(
    f"{'RPM':>8}"
    f"{'Peak flow':>16}"
    f"{'Flow spread':>16}"
    f"{'Ratio spread':>16}"
    f"{'Final bracket':>16}"
    f"{'Residual':>16}"
)

for rpm, rows in sorted(opt_groups.items()):
    rows.sort(key=lambda r: float(r["tolerance"]))

    flows = np.array([
        float(r["peak_flow"]) for r in rows
    ])

    ratios = np.array([
        float(r["peak_ratio"]) for r in rows
    ])

    finest = rows[0]

    print(
        f"{rpm:8.0f}"
        f"{flows[0]:16.9f}"
        f"{np.ptp(flows):16.4e}"
        f"{np.ptp(ratios):16.4e}"
        f"{float(finest['bracket_width']):16.4e}"
        f"{float(finest['residual']):16.4e}"
    )

print("\nSTAGE 7H DERIVATIVE CONVERGENCE")
print("=" * 100)

print(
    f"{'RPM':>8}"
    f"{'h':>12}"
    f"{'First derivative':>20}"
    f"{'Second derivative':>20}"
    f"{'Residual max':>16}"
)

for rpm, rows in sorted(der_groups.items()):
    rows.sort(
        key=lambda r: float(r["h"]),
        reverse=True
    )

    for row in rows:
        residual = max(
            float(row["left_residual"]),
            float(row["center_residual"]),
            float(row["right_residual"])
        )

        print(
            f"{rpm:8.0f}"
            f"{float(row['h']):12.2e}"
            f"{float(row['first_derivative']):20.10e}"
            f"{float(row['second_derivative']):20.10e}"
            f"{residual:16.4e}"
        )

    print("-" * 100)

print("\nSTAGE 7H STABILITY CHECK")
print("=" * 100)

for rpm, rows in sorted(der_groups.items()):
    by_h = {
        float(r["h"]): r
        for r in rows
    }

    # Compare intermediate step sizes, avoiding both
    # large truncation error and extremely small h.
    selected = [2e-3, 1e-3, 5e-4]

    curvature = np.array([
        float(by_h[h]["second_derivative"])
        for h in selected
    ])

    relative_spread = (
        np.ptp(curvature) /
        max(abs(np.mean(curvature)), 1e-30)
    )

    print(
        f"RPM={rpm:.0f} "
        f"curvature_mean={np.mean(curvature):.10f} "
        f"relative_spread={relative_spread:.3e}"
    )
