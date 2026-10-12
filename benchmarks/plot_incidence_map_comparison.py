#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

RESULTS = (
    ROOT / "benchmarks"
    / "results-incidence-map-comparison"
)

DATA = RESULTS / "comparison.csv"

df = pd.read_csv(DATA)

required = {
    "model",
    "mass_flow",
    "rotational_speed_rpm",
    "blade_angle_deg",
    "baseline_corrected_mass_flow",
    "incidence_corrected_mass_flow",
    "baseline_total_pressure_ratio",
    "incidence_total_pressure_ratio",
    "delta_total_pressure_ratio",
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(
        f"Missing CSV columns: {sorted(missing)}"
    )

if len(df) != 105:
    raise RuntimeError(
        f"Expected 105 rows; found {len(df)}"
    )

# Geometry and inlet conditions used by the benchmark.
D1 = 0.15
N_DESIGN = 12000.0
FLOW_DESIGN = 5.0
PREWHIRL = 20.0

def zero_incidence_flow(speed):
    u = np.pi * D1 * speed / 60.0
    u_design = np.pi * D1 * N_DESIGN / 60.0

    return (
        FLOW_DESIGN *
        (u - PREWHIRL) /
        (u_design - PREWHIRL)
    )

RESULTS.mkdir(parents=True, exist_ok=True)

models = [
    "AirLikeIdeal",
    "MethaneIdeal",
    "MethanePR",
]

for model in models:
    subset = df[df.model == model].copy()

    fig, ax = plt.subplots(figsize=(10, 7))

    locus_x = []
    locus_y = []

    for speed, group in subset.groupby(
        "rotational_speed_rpm"
    ):
        group = group.sort_values(
            "baseline_corrected_mass_flow"
        )

        x = group.baseline_corrected_mass_flow.to_numpy()
        baseline = group.baseline_total_pressure_ratio.to_numpy()
        active = group.incidence_total_pressure_ratio.to_numpy()

        baseline_line, = ax.plot(
            x,
            baseline,
            linestyle="--",
            linewidth=1.4,
            label=f"{speed:.0f} rpm baseline",
        )

        ax.plot(
            x,
            active,
            linestyle="-",
            linewidth=1.8,
            color=baseline_line.get_color(),
            label=f"{speed:.0f} rpm incidence",
        )

        # Inlet conditions are fixed within the speed line.
        # Consequently corrected flow is proportional to flow.
        correction_factor = (
            group.baseline_corrected_mass_flow.iloc[0]
            / group.mass_flow.iloc[0]
        )

        flow_zero = zero_incidence_flow(speed)
        x_zero = correction_factor * flow_zero

        # The zero-incidence point need not coincide
        # with a sampled grid point. Interpolate only
        # to display the pressure ratio at that location.
        if x.min() <= x_zero <= x.max():
            y_zero = np.interp(
                x_zero, x, baseline
            )

            locus_x.append(x_zero)
            locus_y.append(y_zero)

    if locus_x:
        order = np.argsort(locus_x)

        ax.plot(
            np.asarray(locus_x)[order],
            np.asarray(locus_y)[order],
            color="black",
            marker="o",
            linewidth=2.0,
            label="Zero-incidence locus",
        )

    ax.set(
        title=f"{model}: incidence-aware compressor map",
        xlabel="Corrected mass flow (kg/s)",
        ylabel="Total pressure ratio",
    )

    ax.grid(True, alpha=0.3)
    ax.legend(
        fontsize=8,
        ncol=2,
        loc="best",
    )

    fig.tight_layout()

    output = RESULTS / f"{model.lower()}_map.png"

    fig.savefig(output, dpi=180)
    plt.close(fig)

    print("Wrote", output)

# Plot the pressure-ratio penalty for all models.
fig, axes = plt.subplots(
    1, 3,
    figsize=(15, 5),
    sharex=True,
    sharey=True,
)

for ax, model in zip(axes, models):
    subset = df[df.model == model]

    speeds = sorted(
        subset.rotational_speed_rpm.unique()
    )

    for speed in speeds:
        group = subset[
            subset.rotational_speed_rpm == speed
        ].sort_values("mass_flow")

        ax.plot(
            group.mass_flow,
            group.delta_total_pressure_ratio,
            marker="o",
            label=f"{speed:.0f} rpm",
        )

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_title(model)
    ax.set_xlabel("Mass flow (kg/s)")
    ax.grid(True, alpha=0.3)

axes[0].set_ylabel(
    "Incidence-induced change in total pressure ratio"
)

axes[-1].legend(fontsize=8)

fig.tight_layout()

output = RESULTS / "incidence_pressure_ratio_penalty.png"

fig.savefig(output, dpi=180)
plt.close(fig)

print("Wrote", output)

print("PASS: incidence compressor-map plots")
