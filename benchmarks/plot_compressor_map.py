#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

root = Path(__file__).resolve().parents[1]
directory = root / "benchmarks/results-compressor-map"

df = pd.read_csv(directory / "map.csv")

for model, group in df.groupby("model"):

    # Plot only continuation results to avoid duplicate points.
    group = group[group["method"] == "continuation"]

    fig, ax = plt.subplots(figsize=(9, 6))

    for speed, line in group.groupby("rotational_speed_rpm"):

        line = line.sort_values("corrected_mass_flow")

        ax.plot(
            line["corrected_mass_flow"],
            line["total_pressure_ratio"],
            marker="o",
            linewidth=1.5,
            label=f"{speed:g} rpm",
        )

    ax.set_xlabel("Corrected mass flow [kg/s]")
    ax.set_ylabel("Total-to-total pressure ratio [-]")

    ax.set_title(
        f"{model}: Compressor Performance Map"
    )

    ax.grid(True, alpha=0.3)
    ax.legend(title="Physical shaft speed")

    fig.tight_layout()

    output = directory / f"{model}_total_pressure_map.png"

    fig.savefig(output, dpi=180)
    plt.close(fig)

    print(f"Generated: {output}")
