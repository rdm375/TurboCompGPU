#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

root = Path(__file__).resolve().parents[1]
directory = root / "benchmarks/results-compressor-map"

df = pd.read_csv(directory / "map.csv")
df = df[df["method"] == "continuation"]

fig, ax = plt.subplots(figsize=(10, 7))

for (model, speed), group in df.groupby(
    ["model", "rotational_speed_rpm"]
):
    group = group.sort_values("corrected_mass_flow")

    ax.plot(
        group["corrected_mass_flow"],
        group["total_pressure_ratio"],
        marker="o",
        label=f"{model}, {speed:g} rpm",
    )

ax.set_xlabel("Corrected mass flow [kg/s]")
ax.set_ylabel("Total-to-total pressure ratio [-]")
ax.set_title("IdealGas vs PengRobinson: Compressor Maps")

ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, ncol=2)

fig.tight_layout()

output = directory / "thermodynamic_model_comparison.png"

fig.savefig(output, dpi=180)
plt.close(fig)

print(f"Generated: {output}")
