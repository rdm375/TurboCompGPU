#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

root = (
    Path(__file__).resolve().parent
    / "results-aerodynamic-transitions"
)


def load(name):
    with (root / name).open(newline="") as f:
        return list(csv.DictReader(f))


transitions = load("transitions.csv")
sensitivities = load("peak_sensitivities.csv")
boundaries = load("boundary_derivatives.csv")

print("\nSTAGE 7I.2 TRANSITION ROOTS")
print("=" * 100)

print(
    f"{'RPM':>7}"
    f"{'Kinc':>9}"
    f"{'Dcrit':>9}"
    f"{'Boundary':>12}"
    f"{'Kdiff root':>16}"
    f"{'Derivative':>18}"
    f"{'Bracket':>16}"
)

for row in sorted(
    transitions,
    key=lambda r: (
        float(r["rpm"]),
        float(r["k_inc"]),
        float(r["d_crit"]),
        r["boundary"]
    )
):
    print(
        f"{float(row['rpm']):7.0f}"
        f"{float(row['k_inc']):9.2f}"
        f"{float(row['d_crit']):9.2f}"
        f"{row['boundary']:>12}"
        f"{float(row['transition_coefficient']):16.9f}"
        f"{float(row['derivative']):18.6e}"
        f"{float(row['bracket_width']):16.3e}"
    )

print("\nTRANSITION COUNTS")
print("=" * 100)

counts = Counter(
    r["boundary"] for r in transitions
)

print("Low-flow crossings:", counts["low"])
print("Capacity crossings:", counts["capacity"])
print("Total crossings:", len(transitions))

print("\nBASELINE PEAK SENSITIVITIES")
print("=" * 100)

print(
    f"{'RPM':>7}"
    f"{'Peak flow':>16}"
    f"{'Slope':>18}"
    f"{'Curvature':>18}"
    f"{'Mixed derivative':>20}"
    f"{'dPeak/dKdiff':>18}"
)

for row in sorted(
    sensitivities,
    key=lambda r: float(r["rpm"])
):
    print(
        f"{float(row['rpm']):7.0f}"
        f"{float(row['peak_flow']):16.9f}"
        f"{float(row['first_derivative']):18.6e}"
        f"{float(row['curvature']):18.6e}"
        f"{float(row['mixed_derivative']):20.6e}"
        f"{float(row['dpeak_dkdiff']):18.9f}"
    )

print("\nBOUNDARY DERIVATIVE SIGN AUDIT")
print("=" * 100)

for rpm in sorted({
    float(r["rpm"]) for r in boundaries
}):
    subset = [
        r for r in boundaries
        if float(r["rpm"]) == rpm
        and float(r["k_inc"]) == 1.0
        and float(r["d_crit"]) == 1.5
    ]

    print(f"\nRPM={rpm:.0f}")

    for kd in sorted({
        float(r["k_diff"]) for r in subset
    }):
        low = next(
            float(r["derivative"])
            for r in subset
            if float(r["k_diff"]) == kd
            and r["boundary"] == "low"
        )

        high = next(
            float(r["derivative"])
            for r in subset
            if float(r["k_diff"]) == kd
            and r["boundary"] == "capacity"
        )

        print(
            f"Kdiff={kd:4.2f} "
            f"low_slope={low:+.8e} "
            f"capacity_slope={high:+.8e}"
        )
