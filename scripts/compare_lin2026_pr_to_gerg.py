#!/usr/bin/env python3
"""Compare the Lin-composition PR property grid with NIST GERG-2008.

Requires the `pyaga8` bindings to the NIST AGA8 GERG-2008 implementation.
The PR grid must include enthalpy_relative_to_reference_J_kg and cp_over_cv,
as written by turbocompgpu_natural_gas_property_grid.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

try:
    from pyaga8 import Composition, Gerg2008
except ImportError as error:
    raise SystemExit(
        "pyaga8 is required (NIST AGA8 GERG-2008 bindings); "
        "install it in the active Python environment first"
    ) from error


COMPOSITION = {
    "methane": 0.9254,
    "nitrogen": 0.0084,
    "carbon_dioxide": 0.0189,
    "ethane": 0.0396,
    "propane": 0.0034,
    "n_butane": 0.0009,
    "isobutane": 0.0012,
    "isopentane": 0.0022,
}
REFERENCE_T_K = 298.15
REFERENCE_P_PA = 101325.0
PAPER_MAX_DEVIATION_PERCENT = {
    "Z": 1.09,
    "density_kg_m3": 1.11,
    "enthalpy_relative_to_reference_J_kg": 1.35,
    "cp_over_cv": 1.72,
}
DIAGNOSTIC_PROPERTIES = ("cp_J_kg_K", "cv_J_kg_K")


def gerg_state(model: Gerg2008, temperature_k: float, pressure_pa: float) -> dict[str, float]:
    model.temperature = temperature_k
    model.pressure = pressure_pa / 1000.0  # pyaga8 takes kPa
    model.calc_density(0)
    model.calc_properties()
    return {
        "Z": float(model.z),
        # GERG density is mol/L and molar mass is g/mol; their product is kg/m3.
        "density_kg_m3": float(model.d * model.mm),
        "enthalpy_J_kg": float(model.h * 1000.0 / model.mm),
        "cp_J_kg_K": float(model.cp * 1000.0 / model.mm),
        "cv_J_kg_K": float(model.cv * 1000.0 / model.mm),
        "cp_over_cv": float(model.cp / model.cv),
        "isentropic_exponent": float(model.kappa),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pr_grid", type=Path, help="PR property-grid CSV")
    parser.add_argument("--output", type=Path, required=True, help="pointwise comparison CSV")
    args = parser.parse_args()

    composition = Composition()
    for component, fraction in COMPOSITION.items():
        setattr(composition, component, fraction)
    gerg = Gerg2008()
    gerg.set_composition(composition)
    reference = gerg_state(gerg, REFERENCE_T_K, REFERENCE_P_PA)

    with args.pr_grid.open(newline="", encoding="utf-8") as stream:
        pr_rows = list(csv.DictReader(stream))
    required = {
        "temperature_K", "pressure_MPa", "Z", "density_kg_m3",
        "enthalpy_relative_to_reference_J_kg", "cp_J_kg_K", "cv_J_kg_K",
        "cp_over_cv",
    }
    if not pr_rows or not required.issubset(pr_rows[0]):
        missing = sorted(required - set(pr_rows[0] if pr_rows else []))
        raise SystemExit(f"PR grid is empty or lacks required columns: {missing}")

    output_rows: list[dict[str, str]] = []
    error_samples: dict[str, list[float]] = {
        key: [] for key in (*PAPER_MAX_DEVIATION_PERCENT, *DIAGNOSTIC_PROPERTIES)
    }
    for row in pr_rows:
        temperature = float(row["temperature_K"])
        pressure_pa = float(row["pressure_MPa"]) * 1.0e6
        ref_state = gerg_state(gerg, temperature, pressure_pa)
        ref_state["enthalpy_relative_to_reference_J_kg"] = (
            ref_state["enthalpy_J_kg"] - reference["enthalpy_J_kg"]
        )
        pr_values = {
            "Z": float(row["Z"]),
            "density_kg_m3": float(row["density_kg_m3"]),
            "enthalpy_relative_to_reference_J_kg": float(
                row["enthalpy_relative_to_reference_J_kg"]
            ),
            "cp_J_kg_K": float(row["cp_J_kg_K"]),
            "cv_J_kg_K": float(row["cv_J_kg_K"]),
            "cp_over_cv": float(row["cp_over_cv"]),
        }
        for key in (*PAPER_MAX_DEVIATION_PERCENT, *DIAGNOSTIC_PROPERTIES):
            reference_value = ref_state[key]
            relative_error = 100.0 * (pr_values[key] - reference_value) / abs(reference_value)
            error_samples[key].append(relative_error)

        comparison = {
            "temperature_K": f"{temperature:.8g}",
            "pressure_MPa": f"{pressure_pa / 1e6:.8g}",
        }
        for key in ("Z", "density_kg_m3", "cp_J_kg_K", "cv_J_kg_K", "cp_over_cv"):
            comparison[f"pr_{key}"] = f"{pr_values[key]:.12g}"
            comparison[f"gerg_{key}"] = f"{ref_state[key]:.12g}"
        comparison["gerg_isentropic_exponent"] = f"{ref_state['isentropic_exponent']:.12g}"
        comparison["pr_enthalpy_relative_J_kg"] = f"{pr_values['enthalpy_relative_to_reference_J_kg']:.12g}"
        comparison["gerg_enthalpy_relative_J_kg"] = f"{ref_state['enthalpy_J_kg'] - reference['enthalpy_J_kg']:.12g}"
        for key in (*PAPER_MAX_DEVIATION_PERCENT, *DIAGNOSTIC_PROPERTIES):
            comparison[f"{key}_relative_error_percent"] = f"{error_samples[key][-1]:.12g}"
        output_rows.append(comparison)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    print("Lin-composition PR vs NIST GERG-2008; 20-point published T,p grid")
    print("Enthalpy is shifted for each EOS to its own value at 298.15 K, 101.325 kPa.")
    print("Specific heat ratio means cp/cv; GERG isentropic exponent is reported separately.")
    print("quantity,max_abs_error_pct,mean_abs_error_pct,paper_max_pct,result")
    for key in (*PAPER_MAX_DEVIATION_PERCENT, *DIAGNOSTIC_PROPERTIES):
        values = error_samples[key]
        max_abs = max(abs(value) for value in values)
        mean_abs = sum(abs(value) for value in values) / len(values)
        if key in PAPER_MAX_DEVIATION_PERCENT:
            limit = PAPER_MAX_DEVIATION_PERCENT[key]
            status = "PASS" if max_abs <= limit else "FAIL"
            limit_text = f"{limit:.2f}"
        else:
            status = "diagnostic"
            limit_text = "not reported"
        print(f"{key},{max_abs:.4f},{mean_abs:.4f},{limit_text},{status}")
    print(f"Pointwise CSV: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
