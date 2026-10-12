#!/usr/bin/env python3
"""Summarize TurboCompGPU diffusion activation-boundary validation."""

import csv
import math
from collections import Counter
from pathlib import Path

SOURCE = Path("benchmarks/results-diffusion-boundary/boundary.csv")
OUTPUT = Path("benchmarks/results-diffusion-boundary/validation_summary.txt")

with SOURCE.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == 15, f"Expected 15 cases, found {len(rows)}"

statuses = Counter(row["status"] for row in rows)
assert statuses["converged"] == 12
assert statuses["not_bracketed"] == 3

converged = [r for r in rows if r["status"] == "converged"]

def values(column):
    result = [float(r[column]) for r in converged]
    assert all(math.isfinite(x) for x in result)
    return result

max_ratio_residual = max(abs(x) for x in values("diffusion_residual"))
max_bracket_width = max(values("bracket_width"))
max_pressure_difference = max(
    abs(x) for x in values("delta_outlet_pressure")
)
max_temperature_difference = max(
    abs(x) for x in values("delta_outlet_temperature")
)
max_coupled_residual = max(
    values("coupled_full_residual_norm")
)

assert max_ratio_residual <= 1e-8
assert max_bracket_width <= 1e-6
assert max_coupled_residual <= 1e-8

lines = [
    "TurboCompGPU M1.6b.3 - Diffusion Activation Boundary",
    "====================================================",
    "",
    "Model:",
    "  Relative-flow diffusion ratio: Dr = W1 / W2",
    "  Critical diffusion ratio: 1.5",
    "  Diffusion coefficient: 1.0 (exploratory)",
    "  Diffusion loss: 0.5*K*W1^2*max(0,Dr-1.5)^2",
    "",
    "Operating conditions:",
    "  Inlet temperature: 300 K",
    "  Inlet pressure: 500000 Pa",
    "  Mass-flow search interval: [2,8] kg/s",
    "  Speeds: 9000,10500,12000,13500,15000 rpm",
    "  Models: AirLikeIdeal, MethaneIdeal, MethanePR",
    "",
    "Validation:",
    f"  Cases: {len(rows)}",
    f"  Converged boundaries: {statuses['converged']}",
    f"  Unbracketed cases: {statuses['not_bracketed']}",
    "  Crossing-direction checks: 12/12 passed",
    "  Coupled nonlinear validation: 12/12 passed",
    f"  Maximum |Dr-1.5|: {max_ratio_residual:.12g}",
    f"  Maximum bracket width (kg/s): {max_bracket_width:.12g}",
    f"  Maximum |delta outlet pressure| (Pa): "
    f"{max_pressure_difference:.12g}",
    f"  Maximum |delta outlet temperature| (K): "
    f"{max_temperature_difference:.12g}",
    f"  Maximum coupled full residual norm: "
    f"{max_coupled_residual:.12g}",
    "",
    "Converged activation boundaries:",
    "  Model            RPM       Mass flow (kg/s)     Phi1",
]

for row in converged:
    lines.append(
        f"  {row['model']:<16} "
        f"{float(row['rotational_speed_rpm']):>7.0f} "
        f"{float(row['critical_mass_flow']):>19.9f} "
        f"{float(row['phi1']):>12.9f}"
    )

lines += [
    "",
    "Interpretation:",
    "  The boundary identifies activation of an exploratory",
    "  relative-flow diffusion-loss correlation.",
    "  It is not an experimentally validated surge boundary.",
    "  Unbracketed means no endpoint sign change over [2,8] kg/s;",
    "  it does not prove that no root exists within or outside",
    "  that interval.",
    "",
]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text("\n".join(lines))

print(OUTPUT.read_text())
