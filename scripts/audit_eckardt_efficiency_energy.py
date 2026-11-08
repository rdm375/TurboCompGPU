#!/usr/bin/env python3
"""Pair Eckardt efficiency digitization with pressure ratio by speedline interpolation.

Pressure ratio and efficiency originate from separate datasets. This script
interpolates the measured pressure map only within each measured speedline and
marks efficiency points outside that range as unpaired. It then compares
implied measured outlet temperature/work with each model's prediction.
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
