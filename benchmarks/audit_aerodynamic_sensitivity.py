#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

root = (
    Path(__file__).resolve().parent
    / "results-aerodynamic-sensitivity"
)

with (root / "sensitivity.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))

print("\nSTAGE 7I CLASSIFICATION AUDIT")
print("=" * 85)

counts = Counter(r["classification"] for r in rows)

print("Total combinations:", len(rows))

for classification, count in sorted(counts.items()):
    print(f"{classification:24s} {count:5d}")

print("\nSTAGE 7I PEAK PREVALENCE BY SPEED")
print("=" * 85)

groups = defaultdict(list)

for row in rows:
    groups[float(row["rpm"])].append(row)

for rpm, group in sorted(groups.items()):
    peaks = [
        r for r in group
        if r["peak_found"] == "1"
    ]

    intervals = np.array([
        float(r["relative_interval"])
        for r in peaks
    ])

    print(
        f"RPM={rpm:7.0f} "
        f"interior={len(peaks):3d}/{len(group):3d} "
        f"interval_min="
        f"{np.min(intervals) if len(intervals) else float('nan'):.6f} "
        f"interval_median="
        f"{np.median(intervals) if len(intervals) else float('nan'):.6f} "
        f"interval_max="
        f"{np.max(intervals) if len(intervals) else float('nan'):.6f}"
    )

print("\nSTAGE 7I BASELINE REPRODUCTION")
print("=" * 85)

reference = {
    9000: 0.517155458,
    10500: 0.817166005,
    12000: 1.140494942,
    13500: 1.494026868,
    15000: 1.884257360,
}

baseline_rows = [
    r for r in rows
    if float(r["k_inc"]) == 1.0
    and float(r["k_diff"]) == 1.0
    and float(r["d_crit"]) == 1.5
]

for row in sorted(
    baseline_rows,
    key=lambda r: float(r["rpm"])
):
    rpm = int(float(row["rpm"]))
    peak = float(row["peak_flow"])
    error = peak - reference[rpm]

    print(
        f"RPM={rpm:5d} "
        f"peak={peak:.9f} "
        f"Stage7H={reference[rpm]:.9f} "
        f"difference={error:+.3e}"
    )

print("\nSTAGE 7I PARAMETER SENSITIVITY")
print("=" * 85)

for parameter in ("k_inc", "k_diff", "d_crit"):
    parameter_groups = defaultdict(list)

    for row in rows:
        parameter_groups[float(row[parameter])].append(row)

    print(f"\n{parameter}")

    for value, group in sorted(parameter_groups.items()):
        peaks = [
            r for r in group
            if r["peak_found"] == "1"
        ]

        intervals = [
            float(r["relative_interval"])
            for r in peaks
        ]

        print(
            f"  {value:5.2f}: "
            f"peaks={len(peaks):3d}/{len(group):3d} "
            f"median_interval="
            f"{np.median(intervals) if intervals else float('nan'):.6f}"
        )
