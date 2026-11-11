#!/usr/bin/env python3
"""Pair Eckardt efficiency digitization with pressure ratio by speedline interpolation.

Pressure ratio and efficiency originate from separate datasets. This script
interpolates the measured pressure map only within each measured speedline and
marks efficiency points outside that range as unpaired. It then compares
implied measured outlet temperature/work with each model's prediction. Select
``--thermodynamic-model pr-air`` to use the benchmark's Peng-Robinson air EOS
for that conversion rather than the ideal-gas diagnostic.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


DEFAULT_PR = Path("data/validation/eckardt_impeller_a.csv")
DEFAULT_EFFICIENCY = Path("data/validation/eckardt_impeller_a_efficiency_digitized.csv")
DEFAULT_MODEL = Path(
    "benchmarks/results-eckardt-impeller-a/impeller_incidence_choke_efficiency.csv"
)
DEFAULT_OUTPUT = Path(
    "benchmarks/results-eckardt-impeller-a/paired_energy_audit.csv"
)


class PengRobinsonAir:
    """Match the pseudo-pure PR air model used by the Eckardt benchmark."""

    def __init__(self, cp: float, molar_mass: float = 0.0289652,
                 critical_temperature: float = 132.53,
                 critical_pressure: float = 3.786e6,
                 acentric_factor: float = 0.035,
                 reference_temperature: float = 298.15,
                 reference_pressure: float = 101325.0):
        self.cp = cp
        self.R = 8.31446261815324 / molar_mass
        self.Tc = critical_temperature
        self.pc = critical_pressure
        self.omega = acentric_factor
        self.Tref = reference_temperature
        self.pref = reference_pressure
        self.kappa = (0.37464 + 1.54226 * self.omega
                      - 0.26992 * self.omega**2)
        self.a0 = 0.45723552892138218938 * self.R**2 * self.Tc**2 / self.pc
        self.b = 0.07779607390388845597 * self.R * self.Tc / self.pc

    def state(self, temperature: float, pressure: float) -> tuple[float, float]:
        """Return mass-specific PR enthalpy and entropy, matching ThermoGPU."""
        if temperature <= 0.0 or pressure <= 0.0:
            raise ValueError("PR state requires positive temperature and pressure")
        root_tr = math.sqrt(temperature / self.Tc)
        alpha_root = 1.0 + self.kappa * (1.0 - root_tr)
        attraction = self.a0 * alpha_root**2
        attraction_derivative = (
            2.0 * self.a0 * alpha_root *
            (-self.kappa / (2.0 * self.Tc * root_tr))
        )

        # The queried states are supercritical vapor states. Bracket the unique
        # physical volume root instead of relying on a cubic-root convention.
        lo = self.b * (1.0 + 1.0e-12)
        hi = max(2.0 * self.R * temperature / pressure, 2.0 * self.b)

        def pressure_at(volume: float) -> float:
            return (self.R * temperature / (volume - self.b)
                    - attraction / (volume * volume + 2.0 * self.b * volume
                                    - self.b * self.b))

        while pressure_at(hi) > pressure:
            hi *= 2.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if pressure_at(mid) > pressure:
                lo = mid
            else:
                hi = mid
        volume = 0.5 * (lo + hi)
        compressibility = pressure * volume / (self.R * temperature)
        sqrt2 = math.sqrt(2.0)
        log_ratio = math.log(
            (volume + (1.0 + sqrt2) * self.b) /
            (volume + (1.0 - sqrt2) * self.b)
        )
        residual_u = (
            (temperature * attraction_derivative - attraction) * log_ratio /
            (2.0 * sqrt2 * self.b)
        )
        residual_s = (
            self.R * math.log(compressibility * (volume - self.b) / volume)
            + attraction_derivative * log_ratio / (2.0 * sqrt2 * self.b)
        )
        enthalpy = (
            self.cp * (temperature - self.Tref) + residual_u
            + pressure * volume - self.R * temperature
        )
        entropy = (
            self.cp * math.log(temperature / self.Tref)
            - self.R * math.log(pressure / self.pref) + residual_s
        )
        return enthalpy, entropy

    def isentropic_enthalpy_rise(
        self, inlet_temperature: float, inlet_pressure: float,
        outlet_pressure: float,
    ) -> float:
        h1, s1 = self.state(inlet_temperature, inlet_pressure)
        lo, hi = 100.0, 1200.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if self.state(mid, outlet_pressure)[1] < s1:
                lo = mid
            else:
                hi = mid
        h_is, _ = self.state(0.5 * (lo + hi), outlet_pressure)
        return h_is - h1

    def temperature_at_enthalpy(
        self, pressure: float, target_enthalpy: float,
    ) -> float:
        lo, hi = 100.0, 1200.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if self.state(mid, pressure)[0] < target_enthalpy:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def read_numeric_csv(path: Path, columns: dict[str, str]) -> list[dict[str, float]]:
    result: list[dict[str, float]] = []
    with path.open(newline="") as stream:
        for row in csv.reader(line for line in stream if not line.lstrip().startswith("#")):
            if not row:
                continue
            # Input validation files are headerless; benchmark result CSVs have headers.
            if row[0].strip().lower() == "model":
                names = [name.strip() for name in row]
                continue
            if "names" in locals():
                values = {name: value for name, value in zip(names, row)}
                result.append({key: float(values[col]) for key, col in columns.items()})
            else:
                vals = [float(value) for value in row]
                result.append({key: vals[int(col)] for key, col in columns.items()})
    return result


def read_model_rows(path: Path) -> list[dict[str, str | float]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(line for line in stream if not line.lstrip().startswith("#"))
        if not reader.fieldnames:
            raise ValueError(f"Missing CSV header in {path}")
        required = {
            "model", "rpm", "mass_flow_kg_s", "measured_eta_tt",
            "predicted_eta_tt", "predicted_total_pressure_ratio",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return [
            {
                "model": row["model"],
                "darcy_factor": float(row.get("darcy_factor") or 0.0),
                "rpm": float(row["rpm"]),
                "mass_flow_kg_s": float(row["mass_flow_kg_s"]),
                "measured_eta_tt": float(row["measured_eta_tt"]),
                "predicted_eta_tt": float(row["predicted_eta_tt"]),
                "predicted_total_pressure_ratio": float(row["predicted_total_pressure_ratio"]),
            }
            for row in reader
        ]


def interpolate_speedline(
    points: list[tuple[float, float]], flow: float
) -> tuple[float | None, float | None, float | None]:
    points = sorted(points)
    for x, y in points:
        if math.isclose(flow, x, rel_tol=0.0, abs_tol=1.0e-10):
            return y, x, x
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 < flow < x1:
            fraction = (flow - x0) / (x1 - x0)
            return y0 + fraction * (y1 - y0), x0, x1
    return None, None, None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pressure-map", type=Path, default=DEFAULT_PR)
    parser.add_argument("--efficiency-map", type=Path, default=DEFAULT_EFFICIENCY)
    parser.add_argument("--model-results", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--inlet-temperature-k", type=float, default=288.15)
    parser.add_argument(
        "--thermodynamic-model", choices=("ideal-gas", "pr-air"),
        default="ideal-gas",
        help=("property model for the work audit; select pr-air for the "
              "Peng-Robinson benchmark runs"),
    )
    parser.add_argument("--cp-j-kg-k", type=float, default=1004.5)
    parser.add_argument("--gas-constant-j-kg-k", type=float, default=287.05)
    args = parser.parse_args()

    pr_rows = read_numeric_csv(
        args.pressure_map,
        {"rpm": "0", "flow": "1", "pr": "2"},
    )
    efficiency_rows = read_numeric_csv(
        args.efficiency_map,
        {"rpm": "0", "flow": "1", "eta": "2"},
    )
    model_rows = read_model_rows(args.model_results)
    if args.inlet_temperature_k <= 0 or args.cp_j_kg_k <= 0:
        raise ValueError("Temperature and cp must be positive")
    gas_gamma = args.cp_j_kg_k / (args.cp_j_kg_k - args.gas_constant_j_kg_k)
    if gas_gamma <= 1.0:
        raise ValueError("cp must be greater than the gas constant")
    pr_air = PengRobinsonAir(args.cp_j_kg_k)

    # Use digitized measured eta as the common operating-point key; model output
    # contains the same digitized samples and may include additional model variants.
    eta_by_point = {(r["rpm"], r["flow"]): r["eta"] for r in efficiency_rows}
    measured_pr: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for row in pr_rows:
        measured_pr[row["rpm"]].append((row["flow"], row["pr"]))

    out_rows: list[dict[str, str | float]] = []
    aggregate: dict[tuple[str, str], list[dict[str, float | bool]]] = defaultdict(list)
    seen_points: set[tuple[float, float]] = set()
    for row in model_rows:
        rpm = float(row["rpm"])
        flow = float(row["mass_flow_kg_s"])
        key = (rpm, flow)
        if key not in eta_by_point:
            continue
        eta_measured = eta_by_point[key]
        pr_measured, flow_lo, flow_hi = interpolate_speedline(measured_pr.get(rpm, []), flow)
        paired = pr_measured is not None
        eta_predicted = float(row["predicted_eta_tt"])
        pr_predicted = float(row["predicted_total_pressure_ratio"])
        darcy_factor = float(row["darcy_factor"])
        row_out: dict[str, str | float] = {
            "model": str(row["model"]),
            "darcy_factor": darcy_factor,
            "rpm": rpm,
            "mass_flow_kg_s": flow,
            "measured_eta_tt_digitized": eta_measured,
            "predicted_eta_tt": eta_predicted,
            "eta_error_percentage_points": 100.0 * (eta_predicted - eta_measured),
            "measured_pressure_ratio_interpolated": "" if pr_measured is None else pr_measured,
            "pressure_ratio_interpolation_flow_low_kg_s": "" if flow_lo is None else flow_lo,
            "pressure_ratio_interpolation_flow_high_kg_s": "" if flow_hi is None else flow_hi,
            "predicted_pressure_ratio": pr_predicted,
            "pressure_ratio_pairing": "linear_interpolation" if paired else "outside_measured_speedline",
        }
        result_stats: dict[str, float | bool] = {"paired": paired}
        if paired:
            assert pr_measured is not None
            if args.thermodynamic_model == "pr-air":
                # Match the benchmark's pseudo-pure dry-air PR parameters,
                # ideal-gas cp, reference state, and 1 atm inlet pressure.
                inlet_pressure = 101325.0
                measured_pressure = inlet_pressure * pr_measured
                predicted_pressure = inlet_pressure * pr_predicted
                h1, _ = pr_air.state(args.inlet_temperature_k, inlet_pressure)
                dh_is_measured = pr_air.isentropic_enthalpy_rise(
                    args.inlet_temperature_k, inlet_pressure, measured_pressure
                )
                dh_is_predicted = pr_air.isentropic_enthalpy_rise(
                    args.inlet_temperature_k, inlet_pressure, predicted_pressure
                )
                dh_measured = dh_is_measured / eta_measured
                dh_predicted = dh_is_predicted / eta_predicted
                t03_measured = pr_air.temperature_at_enthalpy(
                    measured_pressure, h1 + dh_measured
                )
                t03_predicted = pr_air.temperature_at_enthalpy(
                    predicted_pressure, h1 + dh_predicted
                )
            else:
                t03s_measured = args.inlet_temperature_k * pr_measured ** ((gas_gamma - 1.0) / gas_gamma)
                t03_measured = args.inlet_temperature_k + (t03s_measured - args.inlet_temperature_k) / eta_measured
                t03s_predicted = args.inlet_temperature_k * pr_predicted ** ((gas_gamma - 1.0) / gas_gamma)
                t03_predicted = args.inlet_temperature_k + (t03s_predicted - args.inlet_temperature_k) / eta_predicted
                dh_is_measured = args.cp_j_kg_k * (t03s_measured - args.inlet_temperature_k)
                dh_is_predicted = args.cp_j_kg_k * (t03s_predicted - args.inlet_temperature_k)
                dh_measured = args.cp_j_kg_k * (t03_measured - args.inlet_temperature_k)
                dh_predicted = args.cp_j_kg_k * (t03_predicted - args.inlet_temperature_k)
            pr_work_component = (dh_is_predicted - dh_is_measured) / eta_measured
            eta_work_component = dh_is_predicted * (1.0 / eta_predicted - 1.0 / eta_measured)
            row_out.update({
                "implied_measured_T03_K": t03_measured,
                "predicted_T03_K": t03_predicted,
                "T03_residual_predicted_minus_implied_K": t03_predicted - t03_measured,
                "implied_measured_specific_work_J_kg": dh_measured,
                "predicted_specific_work_J_kg": dh_predicted,
                "specific_work_residual_predicted_minus_implied_J_kg": dh_predicted - dh_measured,
                "specific_work_residual_pct_of_implied": 100.0 * (dh_predicted - dh_measured) / dh_measured if dh_measured else math.nan,
                "work_residual_component_from_pressure_ratio_J_kg": pr_work_component,
                "work_residual_component_from_efficiency_J_kg": eta_work_component,
            })
            result_stats.update({
                "pr_error_pct": 100.0 * (pr_predicted - pr_measured) / pr_measured,
                "work_residual": dh_predicted - dh_measured,
                "pr_work_component": pr_work_component,
                "eta_work_component": eta_work_component,
            })
        else:
            row_out.update({
                "implied_measured_T03_K": "",
                "predicted_T03_K": "",
                "T03_residual_predicted_minus_implied_K": "",
                "implied_measured_specific_work_J_kg": "",
                "predicted_specific_work_J_kg": "",
                "specific_work_residual_predicted_minus_implied_J_kg": "",
                "specific_work_residual_pct_of_implied": "",
                "work_residual_component_from_pressure_ratio_J_kg": "",
                "work_residual_component_from_efficiency_J_kg": "",
            })
        out_rows.append(row_out)
        model_name = str(row["model"])
        # Fixed-factor sweeps need separate summaries per factor. The smooth-wall
        # model's factor is calculated per point, so summarize it as one model.
        factor_label = f"{darcy_factor:g}" if model_name == "fixed_galvas_incidence" else "Re-dependent"
        aggregate[(model_name, factor_label)].append(result_stats)
        seen_points.add(key)

    if not out_rows:
        raise ValueError("No model rows matched the digitized efficiency samples")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)

    base_points = len(efficiency_rows)
    print(f"Efficiency digitization points: {base_points}")
    print(f"Energy audit thermodynamics: {args.thermodynamic_model}")
    for (model, factor), rows in sorted(aggregate.items()):
        paired = [r for r in rows if r["paired"]]
        if not paired:
            print(f"{model} (Darcy f={factor}): paired=0/{base_points}")
            continue
        pr_errors = [float(r["pr_error_pct"]) for r in paired]
        work_residuals = [float(r["work_residual"]) for r in paired]
        print(
            f"{model} (Darcy f={factor}): paired={len(paired)}/{base_points}, "
            f"PR MAPE={sum(abs(v) for v in pr_errors)/len(pr_errors):.3f}%, "
            f"work residual mean={sum(work_residuals)/len(work_residuals):.1f} J/kg, "
            f"work residual MAE={sum(abs(v) for v in work_residuals)/len(work_residuals):.1f} J/kg, "
            f"PR component mean={sum(float(r['pr_work_component']) for r in paired)/len(paired):.1f} J/kg, "
            f"eta component mean={sum(float(r['eta_work_component']) for r in paired)/len(paired):.1f} J/kg"
        )
    print(f"Audit CSV: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
