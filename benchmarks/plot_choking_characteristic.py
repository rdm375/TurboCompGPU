#!/usr/bin/env python3
"""Stage 7D: Compressor characteristics with provisional throat capacity."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

LOW_FLOW = (
    ROOT / "benchmarks/results-low-flow-characteristic/"
    "characteristic.csv"
)

ENDPOINTS = (
    ROOT / "benchmarks/results-choking-characteristic/"
    "endpoints.csv"
)

OUT = ROOT / "benchmarks/results-choking-characteristic"
OUT.mkdir(parents=True, exist_ok=True)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


rows = read_csv(LOW_FLOW)
endpoints = read_csv(ENDPOINTS)

assert len(rows) == 105
assert len(endpoints) == 5

by_speed = defaultdict(list)

for row in rows:
    by_speed[float(row["rotational_speed_rpm"])].append(row)

for speed in by_speed:
    by_speed[speed].sort(key=lambda r: float(r["mass_flow"]))

endpoint_by_speed = {
    float(row["rotational_speed_rpm"]): row
    for row in endpoints
}

assert set(by_speed) == set(endpoint_by_speed)


def plot_map(xfield, xlabel, filename):
    fig, ax = plt.subplots(figsize=(11, 8))

    endpoint_x = []
    endpoint_y = []

    for speed, points in sorted(by_speed.items()):
        endpoint = endpoint_by_speed[speed]

        critical_flow = float(endpoint["critical_mass_flow"])

        valid = [
            row for row in points
            if row["mapped"] == "1"
        ]

        below = [
            row for row in valid
            if float(row["mass_flow"]) < critical_flow
        ]

        # Full mathematical operating-point characteristic.
        full_x = [float(r[xfield]) for r in valid]
        full_y = [
            float(r["total_pressure_ratio"])
            for r in valid
        ]

        ax.plot(
            full_x,
            full_y,
            linestyle="--",
            linewidth=1,
            alpha=0.35,
        )

        # Sampled points below the independent capacity.
        x = [float(r[xfield]) for r in below]
        y = [
            float(r["total_pressure_ratio"])
            for r in below
        ]

        # Exact operating-point solution at the predicted
        # choking mass flow, not interpolated map data.
        endpoint_field = (
            "critical_mass_flow"
            if xfield == "mass_flow"
            else xfield
        )
        xe = float(endpoint[endpoint_field])
        ye = float(endpoint["total_pressure_ratio"])

        x.append(xe)
        y.append(ye)

        ax.plot(
            x,
            y,
            marker="o",
            markersize=3,
            linewidth=2,
            label=f"{speed:,.0f} RPM",
        )

        endpoint_x.append(xe)
        endpoint_y.append(ye)

    ax.plot(
        endpoint_x,
        endpoint_y,
        color="black",
        linestyle=":",
        linewidth=2,
        marker="s",
        markersize=7,
        label="Provisional throat-capacity line",
    )

    ax.set(
        xlabel=xlabel,
        ylabel="Total-to-total pressure ratio [-]",
        title="Methane Peng-Robinson compressor characteristic map",
    )

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.text(
        0.5,
        0.01,
        "Endpoint pressure ratios solved directly. "
        "Choking models not yet fully coupled; "
        "surge boundary not modeled.",
        ha="center",
        fontsize=9,
    )

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT / filename, dpi=200)
    plt.close(fig)


plot_map(
    "mass_flow",
    "Mass flow [kg/s]",
    "characteristic_map.png",
)

plot_map(
    "corrected_mass_flow",
    "Corrected mass flow [kg/s]",
    "corrected_characteristic_map.png",
)

print("\nSTAGE 7D ENDPOINT SUMMARY")
print("=" * 100)

print(
    f"{'RPM':>8}"
    f"{'Critical flow':>16}"
    f"{'Pi_tt':>12}"
    f"{'Power kW':>14}"
    f"{'Delta T0':>14}"
    f"{'Delta p0':>14}"
)

for speed, row in sorted(endpoint_by_speed.items()):
    print(
        f"{speed:8.0f}"
        f"{float(row['critical_mass_flow']):16.6f}"
        f"{float(row['total_pressure_ratio']):12.6f}"
        f"{float(row['shaft_power']) / 1000:14.3f}"
        f"{float(row['inlet_total_temperature_difference']):14.6f}"
        f"{float(row['inlet_total_pressure_difference']):14.3f}"
    )

print("\nGenerated:")
print(OUT / "characteristic_map.png")
print(OUT / "corrected_characteristic_map.png")
