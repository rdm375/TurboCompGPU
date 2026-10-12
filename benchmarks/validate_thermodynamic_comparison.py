#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

RESULTS = (
    ROOT
    / "benchmarks"
    / "results-thermodynamic-comparison"
)

DATA = RESULTS / "comparison.csv"
BASELINE = (
    ROOT
    / "benchmarks"
    / "results-compressor-map"
    / "map.csv"
)

df = pd.read_csv(DATA)
baseline = pd.read_csv(BASELINE)

models = {
    "AirLikeIdeal",
    "MethaneIdeal",
    "MethanePR",
}

assert set(df.model) == models

keys = ["mass_flow", "rotational_speed_rpm"]

assert len(df) == 105
assert df.groupby("model").size().eq(35).all()
assert not df.duplicated(["model", *keys]).any()

assert df.converged.eq(1).all()
assert df.validated.eq(1).all()
assert df.basic_physics_valid.eq(1).all()

assert np.isfinite(
    df.select_dtypes(include="number").to_numpy()
).all()

assert np.allclose(
    df.static_pressure_ratio,
    df.outlet_pressure / df.inlet_pressure,
    rtol=1e-10,
    atol=1e-12,
)

assert np.allclose(
    df.total_pressure_ratio,
    df.outlet_total_pressure / df.inlet_total_pressure,
    rtol=1e-10,
    atol=1e-12,
)

assert np.allclose(
    df.shaft_power,
    df.mass_flow * df.specific_work,
    rtol=1e-10,
    atol=1e-8,
)

# M1.6a used these same two thermodynamic models.
# Compare against the independent-solve baseline rows.
comparisons = {
    "AirLikeIdeal": "IdealGas",
    "MethanePR": "PengRobinson",
}

for new_name, old_name in comparisons.items():
    new = (
        df[df.model == new_name]
        .set_index(keys)
        .sort_index()
    )

    old = (
        baseline[
            (baseline.model == old_name)
            & (baseline.method == "independent")
        ]
        .set_index(keys)
        .sort_index()
    )

    assert new.index.equals(old.index)

    for quantity in [
        "outlet_temperature",
        "outlet_pressure",
        "static_pressure_ratio",
        "total_pressure_ratio",
        "specific_work",
        "shaft_power",
        "torque",
        "corrected_mass_flow",
        "corrected_speed_rpm",
    ]:
        assert np.allclose(
            new[quantity],
            old[quantity],
            rtol=1e-8,
            atol=1e-8,
        ), (new_name, quantity)

ideal = (
    df[df.model == "MethaneIdeal"]
    .set_index(keys)
    .sort_index()
)

pr = (
    df[df.model == "MethanePR"]
    .set_index(keys)
    .sort_index()
)

assert ideal.index.equals(pr.index)

# The matched methane models must have the same ideal-gas R.
# At identical T and p, rho_PR / rho_ideal = 1 / Z_PR.
assert np.allclose(
    pr.inlet_density / ideal.inlet_density,
    1.0 / pr.inlet_Z,
    rtol=1e-8,
    atol=1e-10,
)

assert np.allclose(
    pr.outlet_density / ideal.outlet_density
    * ideal.outlet_pressure / pr.outlet_pressure
    * pr.outlet_temperature / ideal.outlet_temperature,
    1.0 / pr.outlet_Z,
    rtol=1e-8,
    atol=1e-10,
)

quantities = [
    "inlet_density",
    "outlet_density",
    "static_pressure_ratio",
    "total_pressure_ratio",
    "specific_work",
]

lines = [
    "M1.6b.1 controlled thermodynamic comparison",
    "=" * 46,
    "105 operating points: PASS",
    "M1.6a baseline comparison: PASS",
    "Pressure-ratio identities: PASS",
    "Power identity: PASS",
    "Matched methane EOS density identities: PASS",
    "",
    "MethanePR minus MethaneIdeal:",
]

for quantity in quantities:
    difference = pr[quantity] - ideal[quantity]

    lines.append(
        f"{quantity}: "
        f"min={difference.min():.10g}, "
        f"max={difference.max():.10g}, "
        f"mean={difference.mean():.10g}"
    )

report = "\n".join(lines) + "\n"

RESULTS.mkdir(parents=True, exist_ok=True)

output = RESULTS / "validation_summary.txt"
output.write_text(report)

print(report)
print("Summary:", output)
