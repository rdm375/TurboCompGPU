#!/usr/bin/env python3
"""Stage 7B low-flow compressor characteristics."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "benchmarks/results-low-flow-characteristic"

with (DIRECTORY / "characteristic.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))

assert len(rows) == 105, f"Expected 105 rows, got {len(rows)}"

groups = defaultdict(list)

for row in rows:
    groups[float(row["rotational_speed_rpm"])].append(row)

for speed in groups:
    groups[speed].sort(key=lambda row: float(row["mass_flow"]))


def value(row, field):
    return float(row[field])


def plot_characteristic(field, ylabel, filename):
    fig, ax = plt.subplots(figsize=(10, 7))

    for speed, points in sorted(groups.items()):
        good = [
            row for row in points
            if row["mapped"] == "1"
            and np.isfinite(value(row, field))
        ]

        flow = [value(row, "mass_flow") for row in good]
        output = [value(row, field) for row in good]

        ax.plot(
            flow,
            output,
            marker="o",
            linewidth=1.7,
            label=f"{speed:,.0f} RPM",
        )

    ax.set(
        xlabel="Mass flow [kg/s]",
        ylabel=ylabel,
        title="Methane PR — low-flow compressor characteristics",
    )

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.text(
        0.5, 0.01,
        "Validated operating-point solutions; "
        "choke admissibility not established",
        ha="center",
        fontsize=9,
    )

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(DIRECTORY / filename, dpi=180)
    plt.close(fig)


plot_characteristic(
    "total_pressure_ratio",
    "Total-to-total pressure ratio [-]",
    "pressure_ratio.png",
)

plot_characteristic(
    "shaft_power",
    "Shaft power [W]",
    "shaft_power.png",
)

fig, ax = plt.subplots(figsize=(10, 7))

for speed, points in sorted(groups.items()):
    flow = np.array([
        value(row, "mass_flow") for row in points
    ])

    margin = np.array([
        value(row, "capacity_margin") for row in points
    ])

    ax.plot(
        flow,
        margin,
        marker="o",
        label=f"{speed:,.0f} RPM",
    )

ax.axhline(0, color="black", linestyle="--")

ax.set(
    xlabel="Mass flow [kg/s]",
    ylabel="Independent throat-capacity margin [-]",
    title="Low-flow choking consistency audit",
)

ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
fig.savefig(DIRECTORY / "capacity_margin.png", dpi=180)
plt.close(fig)

print("\nStage 7B per-speed summary")
print("=" * 76)
print(
    f"{'RPM':>8} {'Mapped':>8} {'Failed':>8} "
    f"{'Below cap.':>12} {'Mapped below':>14}"
)

for speed, points in sorted(groups.items()):
    mapped = sum(row["mapped"] == "1" for row in points)
    below = sum(
        value(row, "capacity_margin") >= 0.0
        for row in points
    )
    mapped_below = sum(
        row["mapped"] == "1"
        and value(row, "capacity_margin") >= 0.0
        for row in points
    )

    print(
        f"{speed:8.0f} {mapped:8d} "
        f"{len(points) - mapped:8d} "
        f"{below:12d} {mapped_below:14d}"
    )

print("\nGenerated plots in:", DIRECTORY)
