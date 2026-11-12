#!/usr/bin/env python3
"""Plot pressure, flow, and configured isothermal temperature vs. x and time."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fields_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("plots"))
    args = parser.parse_args()
    with args.fields_csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("field CSV contains no samples")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    available_legs = {row["leg"] for row in rows}
    legs = [leg for leg in ("upstream", "downstream") if leg in available_legs]
    specs = (
        ("pressure_pa", "Pressure", "Pressure (kPa)", 1.0e-3, "viridis"),
        ("flow_kg_s", "Mass flow", "Mass flow (kg/s)", 1.0, "plasma"),
        ("temperature_k", "Temperature", "Temperature (K)", 1.0, "inferno"),
    )
    for key, title, color_label, scale, palette in specs:
        fig, axes = plt.subplots(
            len(legs), 1, figsize=(10, 4.1 * len(legs)), squeeze=False,
            constrained_layout=True)
        fig.suptitle(f"Two-pipe compressor transient: {title.lower()}")
        for ax, leg in zip(axes[:, 0], legs):
            leg_rows = [row for row in rows if row["leg"] == leg]
            times = sorted({float(row["time_s"]) for row in leg_rows})
            distances = sorted({float(row["distance_m"]) for row in leg_rows})
            ti = {value: i for i, value in enumerate(times)}
            xi = {value: i for i, value in enumerate(distances)}
            values = np.full((len(times), len(distances)), np.nan)
            for row in leg_rows:
                values[ti[float(row["time_s"])], xi[float(row["distance_m"])]] = float(row[key]) * scale
            if np.isnan(values).any():
                raise SystemExit(f"{leg} field samples do not form a complete grid")
            local_distance = np.asarray(distances) - distances[0]
            x, t = np.meshgrid(local_distance, times)
            image = ax.pcolormesh(x, t, values, shading="auto", cmap=palette)
            ax.set_title(f"{leg.capitalize()} pipe")
            ax.set_xlabel("Distance along segment (m)")
            ax.set_ylabel("Time (s)")
            fig.colorbar(image, ax=ax, label=color_label)
            compressor_x = local_distance[-1] if leg == "upstream" else local_distance[0]
            compressor_label = "Compressor inlet" if leg == "upstream" else "Compressor discharge"
            ax.axvline(compressor_x, color="white", linestyle="--", linewidth=1.6)
            ax.annotate(
                compressor_label,
                xy=(compressor_x, 0.97), xycoords=("data", "axes fraction"),
                xytext=(-6 if leg == "upstream" else 6, 0),
                textcoords="offset points",
                ha="right" if leg == "upstream" else "left", va="top",
                color="black", fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 2.0})
        path = args.output_dir / f"pipeline_{key.removesuffix('_pa').removesuffix('_kg_s').removesuffix('_k')}_time_distance.png"
        fig.savefig(path, dpi=170)
        plt.close(fig)
        print(path)


if __name__ == "__main__":
    main()
