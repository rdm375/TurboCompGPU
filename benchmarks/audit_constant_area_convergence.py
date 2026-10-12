#!/usr/bin/env python3

import csv
import math
from pathlib import Path

path = (
    Path(__file__).resolve().parent
    / "results-constant-area-convergence"
    / "constant_area_sensitivity.csv"
)

with path.open(newline="") as f:
    rows = list(csv.DictReader(f))

speeds = [9000, 10500, 12000, 13500, 15000]
steps = [0.04, 0.02, 0.01, 0.005]

def find(rpm, scale):
    matches = [
        r for r in rows
        if float(r["rpm"]) == rpm
        and abs(float(r["scale"]) - scale) < 1e-10
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one row: rpm={rpm}, scale={scale}"
        )

    return matches[0]

def peak(row):
    if row["classification"] != "interior_peak":
        return math.nan

    return float(row["peak_flow"])

print("\nSTAGE 7J.3 ELASTICITY CONVERGENCE")
print("=" * 90)

all_fine_valid = True
all_converged = True

for rpm in speeds:
    baseline = find(rpm, 1.0)
    y0 = peak(baseline)

    print(f"\nRPM={rpm}")
    print(f"{'step':>12}{'elasticity':>20}{'status':>20}")

    estimates = {}

    for h in steps:
        minus = find(rpm, 1.0 - h)
        plus = find(rpm, 1.0 + h)

        ym = peak(minus)
        yp = peak(plus)

        if not all(map(math.isfinite, (y0, ym, yp))):
            print(
                f"{h:12.6f}"
                f"{'nan':>20}"
                f"{'boundary-limited':>20}"
            )
            continue

        elasticity = (yp - ym) / (2.0 * h * y0)
        estimates[h] = elasticity

        print(
            f"{h:12.6f}"
            f"{elasticity:20.10f}"
            f"{'interior':>20}"
        )

    if 0.01 not in estimates or 0.005 not in estimates:
        all_fine_valid = False
        print("Fine-step comparison unavailable")
        continue

    e1 = estimates[0.01]
    e2 = estimates[0.005]

    difference = abs(e2 - e1)

    richardson = (4.0 * e2 - e1) / 3.0

    print(
        "\nFine-step difference:",
        f"{difference:.9e}"
    )

    print(
        "Richardson estimate:",
        f"{richardson:+.10f}"
    )

    # Provisional tolerance, not a solver accuracy bound.
    if difference >= 0.01:
        all_converged = False

print("\nACCEPTANCE")
print("=" * 90)

count_ok = len(rows) == 45

print("Combination count:", len(rows))
print("Expected combinations:", 45)

print("Count:", "PASS" if count_ok else "FAIL")
print(
    "Fine-step availability:",
    "PASS" if all_fine_valid else "REVIEW"
)
print(
    "Fine-step convergence:",
    "PASS" if all_converged else "REVIEW"
)

if not (count_ok and all_fine_valid and all_converged):
    raise SystemExit(2)
