#!/usr/bin/env python3

import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent

def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

constant = read_csv(
    root / "results-constant-area-sensitivity"
    / "constant_area_sensitivity.csv"
)

geometry = read_csv(
    root / "results-geometry-sensitivity"
    / "geometry_sensitivity.csv"
)

speeds = [9000, 10500, 12000, 13500, 15000]

def find(rows, rpm, scale, variable=None):
    matches = [
        r for r in rows
        if float(r["rpm"]) == rpm
        and abs(float(r["scale"]) - scale) < 1e-10
        and (
            variable is None
            or r["variable"] == variable
        )
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one row: {rpm}, {scale}, {variable}"
        )

    return matches[0]

def number(row, field):
    return float(row[field])

def elasticity(plus, minus, baseline, field, delta):
    y0 = number(baseline, field)

    yp = number(plus, field)
    ym = number(minus, field)

    if not all(map(math.isfinite, (y0, yp, ym))):
        return math.nan

    if y0 == 0:
        return math.nan

    return (yp - ym) / (2.0 * delta * y0)

print("\nSTAGE 7J.2 CONSTANT-AREA SENSITIVITY")
print("=" * 105)

print(
    f"{'RPM':>7}"
    f"{'E radius':>13}"
    f"{'E width':>13}"
    f"{'Predicted':>14}"
    f"{'Measured 5%':>15}"
    f"{'Measured 10%':>16}"
    f"{'Difference':>14}"
)

comparisons = []

for rpm in speeds:
    base_g = find(
        geometry, rpm, 1.0, "outlet_radius"
    )

    r_minus = find(
        geometry, rpm, 0.95, "outlet_radius"
    )
    r_plus = find(
        geometry, rpm, 1.05, "outlet_radius"
    )

    b_minus = find(
        geometry, rpm, 0.95, "outlet_width"
    )
    b_plus = find(
        geometry, rpm, 1.05, "outlet_width"
    )

    er = elasticity(
        r_plus, r_minus, base_g, "peak_flow", 0.05
    )

    eb = elasticity(
        b_plus, b_minus, base_g, "peak_flow", 0.05
    )

    predicted = er - eb

    base = find(constant, rpm, 1.0)

    minus5 = find(constant, rpm, 0.95)
    plus5 = find(constant, rpm, 1.05)

    minus10 = find(constant, rpm, 0.90)
    plus10 = find(constant, rpm, 1.10)

    measured5 = elasticity(
        plus5, minus5, base, "peak_flow", 0.05
    )

    measured10 = elasticity(
        plus10, minus10, base, "peak_flow", 0.10
    )

    difference = measured5 - predicted

    comparisons.append({
        "rpm": rpm,
        "predicted": predicted,
        "measured5": measured5,
        "measured10": measured10,
        "difference": difference,
    })

    print(
        f"{rpm:7d}"
        f"{er:13.6f}"
        f"{eb:13.6f}"
        f"{predicted:14.6f}"
        f"{measured5:15.6f}"
        f"{measured10:16.6f}"
        f"{difference:14.6f}"
    )

print("\nGEOMETRY INVARIANTS")
print("=" * 105)

maximum_area_error = max(
    abs(number(r, "relative_area_error"))
    for r in constant
)

print(
    "Maximum relative outlet-area error:",
    f"{maximum_area_error:.9e}"
)

for rpm in speeds:
    baseline = find(constant, rpm, 1.0)
    critical0 = number(baseline, "critical_flow")

    deviations = [
        abs(number(r, "critical_flow") - critical0)
        for r in constant
        if float(r["rpm"]) == rpm
    ]

    print(
        f"RPM={rpm:5d} "
        f"maximum choking-capacity change="
        f"{max(deviations):.9e}"
    )

print("\nCLASSIFICATIONS")
print("=" * 105)

for r in constant:
    if r["classification"] != "interior_peak":
        print(
            f"RPM={r['rpm']} "
            f"scale={r['scale']} "
            f"classification={r['classification']} "
            f"critical_flow={r['critical_flow']}"
        )

print("\nAUDIT SUMMARY")
print("=" * 105)

print("Combinations:", len(constant))

finite_comparisons = [
    r for r in comparisons
    if math.isfinite(r["difference"])
]

if finite_comparisons:
    max_difference = max(
        abs(r["difference"])
        for r in finite_comparisons
    )

    print(
        "Maximum 5% chain-rule discrepancy:",
        f"{max_difference:.6e}"
    )
else:
    print(
        "No finite chain-rule comparisons; "
        "boundary classifications require review."
    )

area_pass = maximum_area_error < 1e-12
count_pass = len(constant) == 25

print(
    "Area invariant:",
    "PASS" if area_pass else "FAIL"
)

print(
    "Combination count:",
    "PASS" if count_pass else "FAIL"
)

if not (area_pass and count_pass):
    raise SystemExit(2)
