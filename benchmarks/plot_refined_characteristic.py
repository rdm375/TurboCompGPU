#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

OUT = ROOT / "benchmarks/results-refined-characteristic"


def read_csv(filename):
    with (OUT / filename).open(newline="") as stream:
        return list(csv.DictReader(stream))


evaluations = read_csv("evaluations.csv")
maxima = read_csv("maxima.csv")

groups = defaultdict(list)

for row in evaluations:
    speed = float(row["rotational_speed_rpm"])
    groups[speed].append(row)

for speed in groups:
    # Optimization traces are not ordered by flow.
    # Duplicate evaluations are harmless for visualization,
    # but removing them makes the plotted curves cleaner.
    unique = {}

    for row in groups[speed]:
        flow = float(row["mass_flow"])
        unique[flow] = row

    groups[speed] = [
        unique[flow]
        for flow in sorted(unique)
    ]

fig, ax = plt.subplots(figsize=(11, 8))

for speed, points in sorted(groups.items()):
    x = [
        float(row["mass_flow"])
        for row in points
    ]

    y = [
        float(row["total_pressure_ratio"])
        for row in points
    ]

    ax.plot(
        x, y,
        linewidth=1.5,
        label=f"{speed:,.0f} RPM"
    )

peak_x = []
peak_y = []

capacity_x = []
capacity_y = []

for row in maxima:
    critical = float(row["critical_mass_flow"])

    capacity_x.append(critical)
    capacity_y.append(
        float(row["critical_total_pressure_ratio"])
    )

    if row["peak_found"] == "1":
        peak_x.append(
            float(row["peak_mass_flow"])
        )

        peak_y.append(
            float(row["peak_total_pressure_ratio"])
        )

ax.plot(
    peak_x,
    peak_y,
    "ko--",
    linewidth=1.5,
    label="Refined characteristic-peak locus"
)

ax.plot(
    capacity_x,
    capacity_y,
    "ks:",
    linewidth=1.5,
    label="Provisional throat-capacity locus"
)

ax.set(
    xlabel="Mass flow [kg/s]",
    ylabel="Total-to-total pressure ratio [-]",
    title=(
        "Methane PR compressor map: "
        "combined incidence and diffusion losses"
    )
)

ax.grid(alpha=0.3)
ax.legend()

fig.text(
    0.5, 0.01,
    "Peak locus is not a validated surge boundary. "
    "Throat-capacity locus is independently predicted.",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=(0, 0.04, 1, 1))

figure = OUT / "refined_operating_envelope.png"
fig.savefig(figure, dpi=200)
plt.close(fig)

print("\nSTAGE 7G REFINED ENVELOPE")
print("=" * 115)

print(
    f"{'RPM':>8}"
    f"{'Peak flow':>14}"
    f"{'Peak Pi_tt':>14}"
    f"{'Capacity':>14}"
    f"{'Capacity Pi_tt':>16}"
    f"{'Flow interval':>16}"
    f"{'Interval %':>14}"
    f"{'Curvature':>16}"
)

for row in maxima:
    speed = float(row["rotational_speed_rpm"])
    peak = float(row["peak_mass_flow"])
    ratio = float(row["peak_total_pressure_ratio"])
    critical = float(row["critical_mass_flow"])
    endpoint_ratio = float(
        row["critical_total_pressure_ratio"]
    )
    interval = float(row["flow_interval"])
    relative = float(row["relative_flow_interval"])
    curvature = float(row["curvature_estimate"])

    print(
        f"{speed:8.0f}"
        f"{peak:14.7f}"
        f"{ratio:14.7f}"
        f"{critical:14.7f}"
        f"{endpoint_ratio:16.7f}"
        f"{interval:16.7f}"
        f"{100 * relative:14.3f}"
        f"{curvature:16.7f}"
    )

print("\nGenerated:", figure)
