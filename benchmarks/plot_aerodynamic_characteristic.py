#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

OUT = ROOT / "benchmarks/results-aerodynamic-characteristic"

with (OUT / "comparison.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))

assert len(rows) == 315, f"Expected 315 rows, got {len(rows)}"

models = ["baseline", "incidence", "combined"]

groups = defaultdict(list)

for row in rows:
    key = (
        row["model"],
        float(row["rotational_speed_rpm"])
    )
    groups[key].append(row)

speeds = sorted({
    float(row["rotational_speed_rpm"])
    for row in rows
})

for key in groups:
    groups[key].sort(key=lambda r: float(r["mass_flow"]))


def valid(row):
    return row["mapped"] == "1"


def number(row, field):
    return float(row[field])


def plot_comparison(field, ylabel, filename):
    fig, axes = plt.subplots(
        len(speeds), 1,
        figsize=(10, 15),
        sharex=True
    )

    for ax, speed in zip(axes, speeds):
        for model in models:
            data = [
                r for r in groups[(model, speed)]
                if valid(r)
            ]

            if not data:
                continue

            x = [number(r, "mass_flow") for r in data]
            y = [number(r, field) for r in data]

            ax.plot(
                x, y,
                marker=".",
                label=model
            )

        ax.set_title(f"{speed:,.0f} RPM")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend()

    axes[-1].set_xlabel("Mass flow [kg/s]")

    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


plot_comparison(
    "total_pressure_ratio",
    "Total pressure ratio",
    "pressure_ratio_comparison.png"
)

plot_comparison(
    "total_entropy_generation",
    "Entropy generation [J/(kg K)]",
    "entropy_comparison.png"
)

plot_comparison(
    "shaft_power",
    "Shaft power [W]",
    "power_comparison.png"
)


print("\nSTAGE 7E CHARACTERISTIC AUDIT")
print("=" * 105)

print(
    f"{'RPM':>8}"
    f"{'Model':>12}"
    f"{'Mapped':>9}"
    f"{'Failed':>9}"
    f"{'Positive slopes':>18}"
    f"{'Negative slopes':>18}"
    f"{'Slope range':>30}"
)

for speed in speeds:
    for model in models:
        data = groups[(model, speed)]
        mapped = [r for r in data if valid(r)]

        slopes = []

        for a, b in zip(mapped[:-1], mapped[1:]):
            dx = number(b, "mass_flow") - number(a, "mass_flow")

            # Do not bridge failed operating points.
            if not np.isclose(dx, 0.1, atol=1e-9):
                continue

            dy = (
                number(b, "total_pressure_ratio") -
                number(a, "total_pressure_ratio")
            )

            slopes.append(dy / dx)

        if slopes:
            slope_range = (
                f"{min(slopes):+.6f} to {max(slopes):+.6f}"
            )
        else:
            slope_range = "N/A"

        print(
            f"{speed:8.0f}"
            f"{model:>12}"
            f"{len(mapped):9d}"
            f"{len(data) - len(mapped):9d}"
            f"{sum(s > 0 for s in slopes):18d}"
            f"{sum(s < 0 for s in slopes):18d}"
            f"{slope_range:>30}"
        )

print("\nGenerated plots in:", OUT)
