#!/usr/bin/env python3

import csv
import math
from collections import Counter
from pathlib import Path

root = (
    Path(__file__).resolve().parent
    / "results-geometry-sensitivity"
)

with (root / "geometry_sensitivity.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))

def get(rpm, variable, scale):
    matches = [
        r for r in rows
        if float(r["rpm"]) == rpm
        and r["variable"] == variable
        and abs(float(r["scale"]) - scale) < 1e-10
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one row: {rpm}, {variable}, {scale}"
        )
    return matches[0]

def value(row, key):
    return float(row[key])

def elasticity(plus, minus, base, key):
    y0 = value(base, key)
    if not math.isfinite(y0) or y0 == 0:
        return math.nan

    return (
        value(plus, key) - value(minus, key)
    ) / (0.10 * y0)

print("\nSTAGE 7J CLASSIFICATION AUDIT")
print("=" * 100)

print("Total combinations:", len(rows))

for name, count in sorted(
    Counter(r["classification"] for r in rows).items()
):
    print(f"{name:25s} {count:3d}")

print("\nSTAGE 7J BASELINE REPRODUCTION")
print("=" * 100)

expected = {
    9000.0: 0.517155458,
    10500.0: 0.817166005,
    12000.0: 1.140494942,
    13500.0: 1.494026868,
    15000.0: 1.884257360,
}

for rpm, reference in expected.items():
    base = get(rpm, "outlet_radius", 1.0)

    peak = value(base, "peak_flow")

    print(
        f"RPM={rpm:5.0f} "
        f"peak={peak:.9f} "
        f"Stage7H={reference:.9f} "
        f"difference={peak-reference:+.3e}"
    )

print("\nSTAGE 7J GEOMETRY ELASTICITIES")
print("=" * 100)

for rpm in expected:
    print(f"\nRPM={rpm:.0f}")

    for variable in ("outlet_radius", "outlet_width"):
        base = get(rpm, "outlet_radius", 1.0)
        minus = get(rpm, variable, 0.95)
        plus = get(rpm, variable, 1.05)

        ep = elasticity(
            plus, minus, base, "peak_flow"
        )

        ec = elasticity(
            plus, minus, base, "critical_flow"
        )

        i0 = value(base, "relative_interval")

        if all(map(math.isfinite, (ep, ec, i0))):
            interval_sensitivity = (1.0 - i0) * (ec - ep)
        else:
            interval_sensitivity = math.nan

        print(
            f"{variable:15s} "
            f"Epeak={ep:+.6f} "
            f"Echoke={ec:+.6f} "
            f"dI/dlnx={interval_sensitivity:+.6f}"
        )

print("\nSTAGE 7J CHOKING CAPACITY VARIATION")
print("=" * 100)

for rpm in expected:
    baseline = get(rpm, "outlet_radius", 1.0)
    critical0 = value(baseline, "critical_flow")

    for variable in ("outlet_radius", "outlet_width"):
        subset = [
            r for r in rows
            if float(r["rpm"]) == rpm
            and r["variable"] == variable
        ]

        deviations = [
            abs(value(r, "critical_flow") - critical0)
            for r in subset
        ]

        print(
            f"RPM={rpm:5.0f} "
            f"{variable:15s} "
            f"max_capacity_change={max(deviations):.6e}"
        )

print("\nSTAGE 7J SOLVER DIAGNOSTICS")
print("=" * 100)

max_choking_residual = max(
    abs(value(r, "choking_mass_residual"))
    for r in rows
)

print(
    "Maximum choking mass residual:",
    f"{max_choking_residual:.6e}"
)

print(
    "Incomplete classifications:",
    sum(
        r["classification"] in (
            "incomplete",
            "refinement_failure",
            "unresolved",
        )
        for r in rows
    )
)
