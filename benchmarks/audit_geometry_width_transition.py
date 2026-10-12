#!/usr/bin/env python3

import csv
import math
from pathlib import Path

path = (
    Path(__file__).resolve().parent
    / "results-geometry-width-transition"
    / "width_transition.csv"
)

with path.open(newline="") as f:
    rows = list(csv.DictReader(f))

rows.sort(key=lambda r: float(r["h"]), reverse=True)

print("\nSTAGE 7J.1 WIDTH TRANSITION CONVERGENCE")
print("=" * 90)

print(
    f"{'h (kg/s)':>12}"
    f"{'width (mm)':>16}"
    f"{'increase (%)':>16}"
    f"{'slope':>16}"
    f"{'bracket (m)':>16}"
)

widths = []

for r in rows:
    h = float(r["h"])
    width = float(r["width_root"])
    increase = 100.0 * float(r["relative_width_change"])
    slope = float(r["derivative"])
    bracket = float(r["bracket_width"])

    widths.append(width)

    print(
        f"{h:12.6f}"
        f"{1000.0 * width:16.9f}"
        f"{increase:16.9f}"
        f"{slope:16.6e}"
        f"{bracket:16.6e}"
    )

print("\nCONVERGENCE DIAGNOSTICS")
print("=" * 90)

spread = max(widths) - min(widths)
fine_spread = max(widths[-3:]) - min(widths[-3:])

print(f"Full width spread: {spread:.9e} m")
print(f"Fine-step spread:  {fine_spread:.9e} m")
print(f"Fine-step spread:  {1000.0*fine_spread:.6e} mm")

for r in rows:
    h = float(r["h"])
    left = float(r["left_derivative"])
    right = float(r["right_derivative"])

    bracketed = (
        math.isfinite(left)
        and math.isfinite(right)
        and left * right < 0
    )

    print(
        f"h={h:.6f} "
        f"left_slope={left:+.6e} "
        f"right_slope={right:+.6e} "
        f"bracketed={bracketed}"
    )

print("\nPROVISIONAL ACCEPTANCE")
print("=" * 90)

# One micrometer in physical outlet width.
tolerance = 1e-6

passed = (
    len(rows) == 5
    and all(math.isfinite(w) for w in widths)
    and fine_spread < tolerance
)

print(
    "Fine-step width convergence:",
    "PASS" if passed else "REVIEW"
)

print(
    "Finest-step transition width:",
    f"{1000.0*widths[-1]:.9f} mm"
)

if not passed:
    raise SystemExit(2)
