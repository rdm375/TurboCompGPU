#!/usr/bin/env python3
"""Plot M1.6b.3 diffusion-aware compressor-map results."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("benchmarks/results-diffusion-map-comparison")
CSV = ROOT / "comparison.csv"
PLOTS = ROOT / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

THRESHOLD = 1.5

MODELS = {
    "AirLikeIdeal": "Air-like ideal gas",
    "MethaneIdeal": "Methane ideal gas",
    "MethanePR": "Methane Peng–Robinson",
}

with CSV.open(newline="") as f:
    rows = list(csv.DictReader(f))

assert len(rows) == 105, f"Expected 105 rows, found {len(rows)}"

data = defaultdict(list)

for row in rows:
    data[row["model"]].append(row)


def value(row, key):
    return float(row[key])


def grouped_by_speed(model):
    groups = defaultdict(list)

    for row in data[model]:
        groups[value(row, "rotational_speed_rpm")].append(row)

    for speed in groups:
        groups[speed].sort(key=lambda r: value(r, "mass_flow"))

    return dict(sorted(groups.items()))


def finish(fig, filename):
    fig.tight_layout()
    path = PLOTS / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(path)


# ------------------------------------------------------------
# 1–3: Relative-flow diffusion ratio
# ------------------------------------------------------------

for model, title in MODELS.items():
    fig, ax = plt.subplots(figsize=(9, 6))

    for speed, points in grouped_by_speed(model).items():
        flow = [value(r, "mass_flow") for r in points]
        ratio = [value(r, "diffusion_ratio") for r in points]

        ax.plot(
            flow,
            ratio,
            marker="o",
            linewidth=1.8,
            label=f"{speed:,.0f} rpm",
        )

    ax.axhline(
        THRESHOLD,
        linestyle="--",
        linewidth=1.5,
        color="black",
        label=f"Diffusion threshold = {THRESHOLD}",
    )

    ax.set(
        title=f"Relative-flow diffusion ratio — {title}",
        xlabel="Mass flow (kg/s)",
        ylabel=r"$D_r = W_1/W_2$",
    )

    ax.grid(True, alpha=0.3)
    ax.legend()

    finish(fig, f"{model}_diffusion_ratio.png")


# ------------------------------------------------------------
# 4–6: Specific diffusion loss
# ------------------------------------------------------------

for model, title in MODELS.items():
    fig, ax = plt.subplots(figsize=(9, 6))

    for speed, points in grouped_by_speed(model).items():
        flow = [value(r, "mass_flow") for r in points]
        loss = [
            value(r, "diffusion_specific_loss")
            for r in points
        ]

        ax.plot(
            flow,
            loss,
            marker="o",
            linewidth=1.8,
            label=f"{speed:,.0f} rpm",
        )

    ax.set(
        title=f"Specific diffusion loss — {title}",
        xlabel="Mass flow (kg/s)",
        ylabel="Specific diffusion loss (J/kg)",
    )

    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    finish(fig, f"{model}_diffusion_loss.png")


# ------------------------------------------------------------
# 7–9: Static compressor pressure ratio
# ------------------------------------------------------------

for model, title in MODELS.items():
    fig, ax = plt.subplots(figsize=(10, 6))

    for speed, points in grouped_by_speed(model).items():
        flow = [value(r, "mass_flow") for r in points]

        incidence = [
            value(r, "incidence_static_pressure_ratio")
            for r in points
        ]

        combined = [
            value(r, "combined_static_pressure_ratio")
            for r in points
        ]

        line, = ax.plot(
            flow,
            incidence,
            linestyle="--",
            linewidth=1.5,
            label=f"{speed:,.0f} rpm — incidence",
        )

        ax.plot(
            flow,
            combined,
            linestyle="-",
            marker="o",
            markersize=4,
            linewidth=1.8,
            color=line.get_color(),
            label=f"{speed:,.0f} rpm — combined",
        )

    ax.set(
        title=f"Compressor pressure ratio — {title}",
        xlabel="Mass flow (kg/s)",
        ylabel="Static pressure ratio",
    )

    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    finish(fig, f"{model}_pressure_ratio.png")


# ------------------------------------------------------------
# 10: Diffusion activation boundary
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

for model, title in MODELS.items():
    speeds = []
    maximum_active_flows = []

    for speed, points in grouped_by_speed(model).items():
        active = [
            value(r, "mass_flow")
            for r in points
            if value(r, "diffusion_specific_loss") > 0
        ]

        if active:
            speeds.append(speed)
            maximum_active_flows.append(max(active))

    ax.plot(
        speeds,
        maximum_active_flows,
        marker="o",
        linewidth=2,
        label=title,
    )

ax.set(
    title="Diffusion activation — sampled operating envelope",
    xlabel="Rotational speed (rpm)",
    ylabel="Maximum sampled mass flow with active diffusion (kg/s)",
)

ax.grid(True, alpha=0.3)
ax.legend()

finish(fig, "diffusion_activation_boundary.png")

print(f"\nPASS: generated 10 plots in {PLOTS}")
