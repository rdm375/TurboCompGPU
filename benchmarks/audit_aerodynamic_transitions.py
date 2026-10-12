#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

root = (
    Path(__file__).resolve().parent
    / "results-aerodynamic-sensitivity"
)

with (root / "sensitivity.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))

for r in rows:
    for key in (
        "rpm", "k_inc", "k_diff", "d_crit",
        "peak_flow", "critical_flow",
        "relative_interval"
    ):
        r[key] = float(r[key])

print("\nCLASSIFICATION BY SPEED AND DIFFUSION COEFFICIENT")
print("=" * 85)

for rpm in sorted({r["rpm"] for r in rows}):
    print(f"\nRPM={rpm:.0f}")

    for kd in sorted({r["k_diff"] for r in rows}):
        subset = [
            r for r in rows
            if r["rpm"] == rpm
            and r["k_diff"] == kd
        ]

        counts = Counter(
            r["classification"] for r in subset
        )

        print(
            f"Kdiff={kd:4.2f} "
            f"low={counts['low_boundary']:2d} "
            f"interior={counts['interior_peak']:2d} "
            f"capacity={counts['capacity_boundary']:2d}"
        )

print("\nDIFFUSION-COEFFICIENT TRANSITIONS")
print("=" * 85)

groups = defaultdict(list)

for r in rows:
    key = (
        r["rpm"],
        r["k_inc"],
        r["d_crit"]
    )

    groups[key].append(r)

transition_counts = Counter()

for key, group in sorted(groups.items()):
    group.sort(key=lambda r: r["k_diff"])

    sequence = tuple(
        r["classification"] for r in group
    )

    transition_counts[sequence] += 1

for sequence, count in transition_counts.most_common():
    print(
        f"{count:3d} cases: "
        + " -> ".join(sequence)
    )

print("\nBASELINE PARAMETER FAMILY")
print("=" * 85)

for rpm in sorted({r["rpm"] for r in rows}):
    print(f"\nRPM={rpm:.0f}")

    subset = [
        r for r in rows
        if r["rpm"] == rpm
        and r["k_inc"] == 1.0
        and r["d_crit"] == 1.5
    ]

    for r in sorted(
        subset,
        key=lambda r: r["k_diff"]
    ):
        print(
            f"Kdiff={r['k_diff']:4.2f} "
            f"class={r['classification']:18s} "
            f"peak={r['peak_flow']:10.6f} "
            f"interval={r['relative_interval']:10.6f}"
        )
