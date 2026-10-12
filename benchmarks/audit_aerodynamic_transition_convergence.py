#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

root = (
    Path(__file__).resolve().parent
    / "results-aerodynamic-transition-convergence"
)


def load(name):
    with (root / name).open(newline="") as f:
        return list(csv.DictReader(f))


roots = load("root_convergence.csv")
sensitivities = load("sensitivity_convergence.csv")

print("\nSTAGE 7I.3 ROOT CONVERGENCE")
print("=" * 100)

groups = defaultdict(list)

for row in roots:
    key = (
        float(row["rpm"]),
        float(row["k_inc"]),
        float(row["d_crit"]),
        row["boundary"],
    )
    groups[key].append(row)

print("Unique transitions:", len(groups))
print("Root evaluations:", len(roots))

max_spread = 0.0

for key, group in sorted(groups.items()):
    values = np.array([
        float(r["k_diff_root"])
        for r in group
    ])

    spread = float(np.ptp(values))
    max_spread = max(max_spread, spread)

    print(
        f"RPM={key[0]:5.0f} "
        f"Kinc={key[1]:4.2f} "
        f"Dcrit={key[2]:4.2f} "
        f"boundary={key[3]:8s} "
        f"root_mean={np.mean(values):.9f} "
        f"spread={spread:.3e}"
    )

print("\nMaximum root spread:", f"{max_spread:.6e}")

print("\nSTAGE 7I.3 PEAK SENSITIVITY CONVERGENCE")
print("=" * 100)

speed_groups = defaultdict(list)

for row in sensitivities:
    speed_groups[float(row["rpm"])].append(row)

for rpm, group in sorted(speed_groups.items()):
    print(f"\nRPM={rpm:.0f}")

    for h in sorted({
        float(r["h"]) for r in group
    }, reverse=True):
        subset = [
            r for r in group
            if float(r["h"]) == h
        ]

        subset.sort(
            key=lambda r: float(r["dk"]),
            reverse=True
        )

        values = np.array([
            float(r["dpeak_dkdiff"])
            for r in subset
        ])

        print(
            f"h={h:.4f} "
            f"sensitivity_mean={np.mean(values):.9f} "
            f"dk_spread={np.ptp(values):.3e}"
        )

    central = next(
        r for r in group
        if float(r["h"]) == 0.001
        and float(r["dk"]) == 0.005
    )

    print(
        "Reference sensitivity:",
        f"{float(central['dpeak_dkdiff']):.9f}"
    )

print("\nSTAGE 7I.3 NUMERICAL ACCEPTANCE")
print("=" * 100)

root_pass = max_spread < 1e-4

print(
    "Root convergence:",
    "PASS" if root_pass else "REVIEW"
)

for rpm, group in sorted(speed_groups.items()):
    selected = [
        float(r["dpeak_dkdiff"])
        for r in group
        if float(r["h"]) in (0.002, 0.001, 0.0005)
        and float(r["dk"]) in (0.01, 0.005, 0.0025)
    ]

    spread = np.ptp(selected)
    mean = np.mean(selected)

    relative_spread = (
        spread / abs(mean)
        if mean != 0 else float("inf")
    )

    print(
        f"RPM={rpm:.0f} "
        f"relative_sensitivity_spread="
        f"{relative_spread:.3e} "
        f"status="
        f"{'PASS' if relative_spread < 1e-3 else 'REVIEW'}"
    )
