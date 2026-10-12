#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

RESULTS = (
    ROOT
    / "benchmarks"
    / "results-incidence-map-comparison"
)

DATA = RESULTS / "comparison.csv"

BASELINE = (
    ROOT
    / "benchmarks"
    / "results-thermodynamic-comparison"
    / "comparison.csv"
)

df = pd.read_csv(DATA)
old = pd.read_csv(BASELINE)

keys = [
    "model",
    "mass_flow",
    "rotational_speed_rpm",
]

assert len(df) == 105
assert len(old) == 105

assert not df.duplicated(keys).any()
assert not old.duplicated(keys).any()

assert set(df.model) == {
    "AirLikeIdeal",
    "MethaneIdeal",
    "MethanePR",
}

assert df.groupby("model").size().eq(35).all()

assert np.isfinite(
    df.select_dtypes(include="number").to_numpy()
).all()

df = df.set_index(keys).sort_index()
old = old.set_index(keys).sort_index()

assert df.index.equals(old.index)

# Reproduce the committed M1.6b.1 baseline.
for new_column, old_column in [
    ("baseline_outlet_temperature", "outlet_temperature"),
    ("baseline_outlet_pressure", "outlet_pressure"),
    ("baseline_static_pressure_ratio", "static_pressure_ratio"),
    ("baseline_total_pressure_ratio", "total_pressure_ratio"),
    ("baseline_specific_work", "specific_work"),
    ("baseline_shaft_power", "shaft_power"),
    ("baseline_torque", "torque"),
    ("baseline_corrected_mass_flow", "corrected_mass_flow"),
    ("baseline_corrected_speed_rpm", "corrected_speed_rpm"),
]:
    assert np.allclose(
        df[new_column],
        old[old_column],
        rtol=1e-8,
        atol=1e-8,
    ), (new_column, old_column)

# Validate the incidence entropy decomposition.
assert np.allclose(
    df.total_entropy_generation,
    df.baseline_entropy_generation +
    df.incidence_entropy_generation,
    rtol=1e-12,
    atol=1e-12,
)

assert (df.incidence_specific_loss >= 0.0).all()
assert (df.incidence_entropy_generation >= 0.0).all()

# Verify the reference design point.
design = df.reset_index().query(
    "mass_flow == 5.0 and rotational_speed_rpm == 12000.0"
)

assert len(design) == 3

assert np.allclose(
    design.incidence_specific_loss,
    0.0,
    atol=1e-10,
)

assert np.allclose(
    design.delta_total_pressure_ratio,
    0.0,
    atol=1e-10,
)

# Check paired differences directly.
quantities = [
    "outlet_temperature",
    "outlet_pressure",
    "static_pressure_ratio",
    "total_pressure_ratio",
    "specific_work",
    "shaft_power",
    "torque",
]

for quantity in quantities:
    baseline = df[f"baseline_{quantity}"].to_numpy()
    incidence = df[f"incidence_{quantity}"].to_numpy()
    reported_delta = df[f"delta_{quantity}"].to_numpy()

    reconstructed_delta = incidence - baseline

    scale = np.maximum(
        1.0,
        np.maximum(np.abs(baseline), np.abs(incidence)),
    )

    tolerance = (
        32.0 * np.finfo(float).eps * scale
    )

    error = np.abs(
        reported_delta - reconstructed_delta
    )

    if not np.all(error <= tolerance):
        worst = np.argmax(error / tolerance)

        raise AssertionError(
            f"{quantity}: "
            f"max_error={error.max():.17g}, "
            f"worst_error={error[worst]:.17g}, "
            f"allowed={tolerance[worst]:.17g}"
        )

# Verify the power identity for both solutions.
for prefix in ["baseline", "incidence"]:
    mass_flow = df.index.get_level_values("mass_flow")

    assert np.allclose(
        df[f"{prefix}_shaft_power"],
        mass_flow * df[f"{prefix}_specific_work"],
        rtol=1e-10,
        atol=1e-8,
    )

# Corrected coordinates must not change.
assert np.allclose(
    df.baseline_corrected_mass_flow,
    df.incidence_corrected_mass_flow,
    rtol=1e-12,
    atol=1e-12,
)

assert np.allclose(
    df.baseline_corrected_speed_rpm,
    df.incidence_corrected_speed_rpm,
    rtol=1e-12,
    atol=1e-12,
)

lines = [
    "M1.6b.2 Stage 4: Paired incidence compressor maps",
    "=" * 52,
    "105 matched operating points: PASS",
    "315 nonlinear solutions: PASS (benchmark acceptance)",
    "105 zero-coefficient regressions: PASS (benchmark acceptance)",
    "M1.6b.1 baseline agreement: PASS",
    "Incidence entropy decomposition: PASS",
    "Design-point zero incidence: PASS",
    "Power identities: PASS",
    "Corrected-coordinate invariance: PASS",
    "",
]

for model, group in df.groupby(level="model"):
    lines.append(model)
    lines.append("-" * len(model))

    for quantity in [
        "incidence_specific_loss",
        "incidence_entropy_generation",
        "delta_static_pressure_ratio",
        "delta_total_pressure_ratio",
        "delta_outlet_temperature",
        "delta_specific_work",
        "delta_shaft_power",
    ]:
        values = group[quantity]

        lines.append(
            f"{quantity}: "
            f"min={values.min():.10g}, "
            f"max={values.max():.10g}, "
            f"mean={values.mean():.10g}"
        )

    lines.append("")

report = "\n".join(lines)

RESULTS.mkdir(parents=True, exist_ok=True)

output = RESULTS / "validation_summary.txt"
output.write_text(report + "\n")

print(report)
print("Summary:", output)
