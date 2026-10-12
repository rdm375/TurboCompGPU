#!/usr/bin/env python3
"""Stage 7A: compressor characteristics and independent choking audit.

The choking boundary is diagnostic only. The compressor operating-point
model does not yet enforce passage choking or establish surge stability.
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmarks/results-compressor-characteristic"
OUT.mkdir(parents=True, exist_ok=True)

MAP_CSV = ROOT / "benchmarks/results-compressor-map/map.csv"
CHOKE_CSV = ROOT / "benchmarks/results-choking-continuation/continuation.csv"

if not MAP_CSV.exists():
    raise SystemExit(f"Missing map dataset: {MAP_CSV}")

if not CHOKE_CSV.exists():
    raise SystemExit(f"Missing choking dataset: {CHOKE_CSV}")


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


map_rows = read_csv(MAP_CSV)
choke_rows = read_csv(CHOKE_CSV)

# Avoid duplicate operating points from the independent/continuation
# comparison. Use the independently validated continuation results.
map_rows = [
    row for row in map_rows
    if row["model"] == "PengRobinson"
    and row["method"] == "continuation"
    and row["converged"] == "1"
    and row["validated"] == "1"
    and row["basic_physics_valid"] == "1"
]

# Choking results use methane Peng-Robinson with prescribed passage turning.
choke = sorted(
    (
        float(row["parameter"]),
        float(row["mdot_critical"])
    )
    for row in choke_rows
    if row["sweep"] == "rpm"
)

if len(choke) < 2:
    raise SystemExit("Insufficient choking-speed data")

if not map_rows:
    raise SystemExit("No validated Peng-Robinson map points")

choke_speed = np.array([x[0] for x in choke])
choke_flow = np.array([x[1] for x in choke])

by_speed = defaultdict(list)

for row in map_rows:
    speed = float(row["rotational_speed_rpm"])
    by_speed[speed].append(row)

for speed in by_speed:
    by_speed[speed].sort(key=lambda r: float(r["mass_flow"]))

# ------------------------------------------------------------
# Diagnostic CSV
# ------------------------------------------------------------

audit_csv = OUT / "characteristic_audit.csv"

with audit_csv.open("w", newline="") as stream:
    writer = csv.writer(stream)

    writer.writerow([
        "rotational_speed_rpm",
        "mass_flow",
        "total_pressure_ratio",
        "static_pressure_ratio",
        "shaft_power",
        "critical_mass_flow_interpolated",
        "capacity_margin",
        "exceeds_independent_throat_capacity",
        "choke_admissibility_established",
    ])

    for speed, rows in sorted(by_speed.items()):
        if not choke_speed[0] <= speed <= choke_speed[-1]:
            raise SystemExit(
                f"Speed {speed} outside choking interpolation range"
            )

        limit = float(np.interp(speed, choke_speed, choke_flow))

        for row in rows:
            flow = float(row["mass_flow"])
            margin = 1.0 - flow / limit

            writer.writerow([
                speed,
                flow,
                float(row["total_pressure_ratio"]),
                float(row["static_pressure_ratio"]),
                float(row["shaft_power"]),
                limit,
                margin,
                int(margin < 0.0),
                0,
            ])

# ------------------------------------------------------------
# Pressure-ratio characteristics
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 7))

for speed, rows in sorted(by_speed.items()):
    flow = np.array([float(r["mass_flow"]) for r in rows])
    ratio = np.array([
        float(r["total_pressure_ratio"]) for r in rows
    ])

    ax.plot(
        flow,
        ratio,
        marker="o",
        linewidth=1.8,
        label=f"{speed:,.0f} RPM",
    )

ax.set(
    xlabel="Mass flow [kg/s]",
    ylabel="Total-to-total pressure ratio [-]",
    title="Methane Peng-Robinson compressor characteristics",
)

ax.grid(True, alpha=0.3)
ax.legend(title="Rotational speed")

fig.text(
    0.5, 0.01,
    "Operating-point solutions; choke admissibility not established",
    ha="center",
    fontsize=9,
)

fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(OUT / "pressure_ratio.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# Power characteristics
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 7))

for speed, rows in sorted(by_speed.items()):
    flow = np.array([float(r["mass_flow"]) for r in rows])
    power = np.array([
        float(r["shaft_power"]) / 1000.0 for r in rows
    ])

    ax.plot(
        flow,
        power,
        marker="o",
        linewidth=1.8,
        label=f"{speed:,.0f} RPM",
    )

ax.set(
    xlabel="Mass flow [kg/s]",
    ylabel="Shaft power [kW]",
    title="Calculated compressor shaft power",
)

ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
fig.savefig(OUT / "shaft_power.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# Independent choking capacity
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(
    choke_speed,
    choke_flow,
    marker="o",
    linewidth=2,
    label="Stage 6H choking capacity",
)

ax.set(
    xlabel="Rotational speed [RPM]",
    ylabel="Critical mass flow [kg/s]",
    title="Independent passage-choking boundary",
)

ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
fig.savefig(OUT / "choking_boundary.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# Capacity-margin audit
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 7))

for speed, rows in sorted(by_speed.items()):
    flow = np.array([float(r["mass_flow"]) for r in rows])
    limit = float(np.interp(speed, choke_speed, choke_flow))

    margin = 1.0 - flow / limit

    ax.plot(
        flow,
        margin,
        marker="o",
        linewidth=1.8,
        label=f"{speed:,.0f} RPM",
    )

ax.axhline(0.0, color="black", linestyle="--")

ax.set(
    xlabel="Mass flow [kg/s]",
    ylabel="Independent throat-capacity margin [-]",
    title="Operating-point / choking-model consistency audit",
)

ax.grid(True, alpha=0.3)
ax.legend()

fig.text(
    0.5, 0.01,
    "Negative margin indicates a model discrepancy, not validated "
    "compressor inadmissibility",
    ha="center",
    fontsize=9,
)

fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(OUT / "capacity_margin.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# Numerical summary
# ------------------------------------------------------------

print("\nStage 7A: compressor characteristic audit")
print("=" * 76)

print(
    f"{'RPM':>8} {'Choke kg/s':>14} "
    f"{'Map min':>12} {'Map max':>12} {'Overlap':>12}"
)

for speed, rows in sorted(by_speed.items()):
    limit = float(np.interp(speed, choke_speed, choke_flow))

    flows = [float(r["mass_flow"]) for r in rows]

    overlap = min(max(flows), limit) - min(flows)

    print(
        f"{speed:8.0f} {limit:14.6f} "
        f"{min(flows):12.6f} {max(flows):12.6f} "
        f"{max(0.0, overlap):12.6f}"
    )

print("\nGenerated:")
for path in sorted(OUT.iterdir()):
    print(" ", path)

print("\nStage 7A diagnostic generation complete.")
print("Choking admissibility remains unresolved.")
