#!/usr/bin/env python3

import csv
import math
from pathlib import Path

path = (
    Path(__file__).resolve().parent
    / "results-constant-area-transition"
    / "constant_area_transition.csv"
)

with path.open(newline="") as f:
    rows = list(csv.DictReader(f))

rows.sort(key=lambda r: float(r["h"]), reverse=True)

print("\nSTAGE 7J.3 CONSTANT-AREA TRANSITION AUDIT")
print("=" * 100)

print(
    f"{'h':>10}"
    f"{'scale':>16}"
    f"{'radius (mm)':>18}"
    f"{'width (mm)':>18}"
    f"{'slope':>18}"
)

scales = []
radii = []
area_errors = []
brackets_valid = True

for row in rows:
    h = float(row["h"])
    scale = float(row["scale_root"])
    radius = float(row["radius_root"])
    width = float(row["width_root"])
    slope = float(row["derivative"])

    left = float(row["left_derivative"])
    right = float(row["right_derivative"])

    scales.append(scale)
    radii.append(radius)

    # Baseline radius = 0.2 m, width = 0.025 m.
    area_error = abs(
        radius * width / (0.2 * 0.025) - 1.0
    )

    area_errors.append(area_error)

    brackets_valid &= (
        math.isfinite(left)
        and math.isfinite(right)
        and left * right < 0.0
    )

    print(
        f"{h:10.6f}"
        f"{scale:16.10f}"
        f"{1000 * radius:18.9f}"
        f"{1000 * width:18.9f}"
        f"{slope:18.6e}"
    )

fine_scale_spread = (
    max(scales[-3:]) - min(scales[-3:])
)

fine_radius_spread = (
    max(radii[-3:]) - min(radii[-3:])
)

max_area_error = max(area_errors)

print("\nCONVERGENCE DIAGNOSTICS")
print("=" * 100)

print(
    "Fine-step scale spread:",
    f"{fine_scale_spread:.9e}"
)

print(
    "Fine-step radius spread:",
    f"{1000 * fine_radius_spread:.9e} mm"
)

print(
    "Maximum relative area error:",
    f"{max_area_error:.9e}"
)

print("All roots bracketed:", brackets_valid)

passed = (
    len(rows) == 5
    and brackets_valid
    and all(math.isfinite(s) for s in scales)
    and fine_scale_spread < 1e-6
    and max_area_error < 1e-12
)

print("\nACCEPTANCE")
print("=" * 100)
print("Result:", "PASS" if passed else "REVIEW")

if not passed:
    raise SystemExit(2)
