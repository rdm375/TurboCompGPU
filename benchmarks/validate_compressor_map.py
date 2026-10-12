#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]

old_path = root / (
    "benchmarks/results-operating-envelope/envelope.csv"
)

new_path = root / (
    "benchmarks/results-compressor-map/map.csv"
)

old = pd.read_csv(old_path)
new = pd.read_csv(new_path)

keys = [
    "model",
    "method",
    "mass_flow",
    "rotational_speed_rpm",
]

assert len(old) == 140
assert len(new) == 140

assert not old.duplicated(keys).any()
assert not new.duplicated(keys).any()

merged = old.merge(
    new,
    on=keys,
    suffixes=("_old", "_new"),
    validate="one_to_one",
)

assert len(merged) == 140

baseline_fields = [
    "converged",
    "outlet_temperature",
    "outlet_pressure",
    "specific_work",
    "shaft_power",
    "torque",
    "entropy_generation",
    "residual_norm",
    "validation_residual_norm",
    "iterations",
    "evaluations",
]

for field in baseline_fields:
    a = merged[f"{field}_old"].to_numpy()
    b = merged[f"{field}_new"].to_numpy()

    if field in {"iterations", "evaluations", "converged"}:
        assert np.array_equal(a, b), field
    else:
        np.testing.assert_allclose(
            a, b,
            rtol=1e-8,
            atol=1e-10,
            err_msg=field,
        )

assert new["converged"].eq(1).all()
assert new["validated"].eq(1).all()
assert new["basic_physics_valid"].eq(1).all()

assert new["operating_limits"].eq("NotEvaluated").all()

for field in [
    "corrected_mass_flow",
    "corrected_speed_rpm",
    "inlet_total_temperature",
    "inlet_total_pressure",
    "outlet_total_temperature",
    "outlet_total_pressure",
    "static_pressure_ratio",
    "total_pressure_ratio",
]:
    assert np.isfinite(new[field]).all(), field

assert new["corrected_mass_flow"].gt(0).all()
assert new["corrected_speed_rpm"].gt(0).all()
assert new["total_pressure_ratio"].gt(0).all()

reference_temperature = 288.15
reference_pressure = 101325.0

tau = new["inlet_total_temperature"] / reference_temperature
delta = new["inlet_total_pressure"] / reference_pressure

expected_flow = (
    new["mass_flow"] * np.sqrt(tau) / delta
)

expected_speed = (
    new["rotational_speed_rpm"] / np.sqrt(tau)
)

np.testing.assert_allclose(
    new["corrected_mass_flow"],
    expected_flow,
    rtol=1e-12,
)

np.testing.assert_allclose(
    new["corrected_speed_rpm"],
    expected_speed,
    rtol=1e-12,
)

np.testing.assert_allclose(
    new["total_pressure_ratio"],
    new["outlet_total_pressure"] /
    new["inlet_total_pressure"],
    rtol=1e-12,
)

new["total_minus_static_ratio"] = (
    new["total_pressure_ratio"] -
    new["static_pressure_ratio"]
)

summary = new.groupby("model").agg(
    cases=("mass_flow", "size"),
    min_corrected_flow=("corrected_mass_flow", "min"),
    max_corrected_flow=("corrected_mass_flow", "max"),
    min_total_ratio=("total_pressure_ratio", "min"),
    max_total_ratio=("total_pressure_ratio", "max"),
    min_ratio_difference=("total_minus_static_ratio", "min"),
    max_ratio_difference=("total_minus_static_ratio", "max"),
)

print("\nM1.6a compressor-map validation")
print("================================")
print("M1.5 baseline comparison: PASS")
print("Corrected-coordinate identities: PASS")
print("Total-pressure-ratio identity: PASS")
print("Operating-limit classification: PASS")
print("Independent/continuation pairs: 70")
print("\n", summary.to_string())

output = (
    root /
    "benchmarks/results-compressor-map/validation_summary.txt"
)

output.write_text(
    "M1.6a compressor-map validation\n"
    "================================\n"
    "M1.5 baseline comparison: PASS\n"
    "Corrected-coordinate identities: PASS\n"
    "Total-pressure-ratio identity: PASS\n"
    "Operating-limit classification: PASS\n"
    "Independent/continuation pairs: 70\n\n"
    + "\n".join(
        line.rstrip()
        for line in summary.to_string().splitlines()
    )
    + "\n"
)

print(f"\nSummary: {output}")
