#!/usr/bin/env python3
"""Summarize the C++ Eckardt one-at-a-time loss ablation CSV."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(line for line in stream if not line.lstrip().startswith("#")))


def interpolate(points: list[tuple[float, float]], flow: float) -> float | None:
    points = sorted(points)
    for x, y in points:
        if math.isclose(x, flow, rel_tol=0.0, abs_tol=1e-10):
            return y
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 < flow < x1:
            w = (flow - x0) / (x1 - x0)
            return y0 + w * (y1 - y0)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ablation-csv", type=Path,
        default=Path("benchmarks/results-eckardt-impeller-a/loss_term_ablation.csv"),
    )
    parser.add_argument(
        "--pressure-map", type=Path,
        default=Path("data/validation/eckardt_impeller_a.csv"),
    )
    parser.add_argument(
        "--efficiency-map", type=Path,
        default=Path("data/validation/eckardt_impeller_a_efficiency_digitized.csv"),
    )
    parser.add_argument("--inlet-temperature-k", type=float, default=288.15)
    parser.add_argument("--cp-j-kg-k", type=float, default=1004.5)
    parser.add_argument("--gas-constant-j-kg-k", type=float, default=287.05)
    args = parser.parse_args()

    pressure: dict[float, list[tuple[float, float]]] = defaultdict(list)
    with args.pressure_map.open(newline="") as stream:
        for fields in csv.reader(line for line in stream if not line.lstrip().startswith("#")):
            if fields:
                pressure[float(fields[0])].append((float(fields[1]), float(fields[2])))

    efficiency: dict[tuple[float, float], float] = {}
    with args.efficiency_map.open(newline="") as stream:
        for fields in csv.reader(line for line in stream if not line.lstrip().startswith("#")):
            if fields:
                efficiency[(float(fields[0]), float(fields[1]))] = float(fields[2])

    gamma = args.cp_j_kg_k / (args.cp_j_kg_k - args.gas_constant_j_kg_k)
    exponent = (gamma - 1.0) / gamma
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(args.ablation_csv):
        grouped[row["scenario"]].append(row)
    baseline_name = (
        "all_losses_enthalpy_equivalent"
        if "all_losses_enthalpy_equivalent" in grouped
        else "all_losses"
    )
    if baseline_name not in grouped:
        raise ValueError("Ablation CSV has no recognized full-loss baseline rows")

    baseline_by_point = {
        (float(row["rpm"]), float(row["mass_flow_kg_s"])): row
        for row in grouped[baseline_name]
    }
    summaries = {}
    print(
        "scenario,eta_bias_pp,eta_MAE_pp,eta_max_abs_pp,"
        "pressure_ratio_MAPE_pct,paired_work_residual_mean_kJ_kg,"
        "paired_work_residual_MAE_kJ_kg,delta_mean_eta_vs_full_pp,"
        "delta_mean_PR_vs_full_pct,delta_mean_work_vs_full_kJ_kg,paired_points"
    )
    for scenario, rows in sorted(grouped.items()):
        eta_errors: list[float] = []
        pr_errors: list[float] = []
        work_residuals: list[float] = []
        delta_eta: list[float] = []
        delta_pr: list[float] = []
        delta_work: list[float] = []
        paired = 0
        for row in rows:
            rpm = float(row["rpm"])
            flow = float(row["mass_flow_kg_s"])
            eta_measured = efficiency[(rpm, flow)]
            eta_predicted = float(row["predicted_eta_tt"])
            pr_predicted = float(row["predicted_total_pressure_ratio"])
            eta_error = 100.0 * (eta_predicted - eta_measured)
            eta_errors.append(eta_error)
            base = baseline_by_point[(rpm, flow)]
            delta_eta.append(eta_error - 100.0 * (
                float(base["predicted_eta_tt"]) - eta_measured))
            base_pr = float(base["predicted_total_pressure_ratio"])
            delta_pr.append(100.0 * (pr_predicted - base_pr) / base_pr)
            delta_work.append(
                float(row["impeller_specific_work_j_kg"])
                - float(base["impeller_specific_work_j_kg"])
            )
            pr_measured = interpolate(pressure.get(rpm, []), flow)
            if pr_measured is None:
                continue
            paired += 1
            pr_errors.append(100.0 * (pr_predicted - pr_measured) / pr_measured)
            t03s_measured = args.inlet_temperature_k * pr_measured**exponent
            t03_measured = args.inlet_temperature_k + (
                t03s_measured - args.inlet_temperature_k
            ) / eta_measured
            t03_model = float(row["outlet_total_temperature_k"])
            work_residuals.append(args.cp_j_kg_k * (t03_model - t03_measured))

        def mean(values: list[float]) -> float:
            return sum(values) / len(values) if values else math.nan

        summary = {
            "eta_bias": mean(eta_errors),
            "eta_mae": mean([abs(x) for x in eta_errors]),
            "eta_max": max(abs(x) for x in eta_errors),
            "pr_mape": mean([abs(x) for x in pr_errors]),
            "work_bias": mean(work_residuals) / 1000.0,
            "work_mae": mean([abs(x) for x in work_residuals]) / 1000.0,
            "delta_eta": mean(delta_eta),
            "delta_pr": mean(delta_pr),
            "delta_work": mean(delta_work) / 1000.0,
            "paired": paired,
        }
        summaries[scenario] = summary
        print(
            f"{scenario},{summary['eta_bias']:.3f},{summary['eta_mae']:.3f},"
            f"{summary['eta_max']:.3f},{summary['pr_mape']:.3f},"
            f"{summary['work_bias']:.3f},{summary['work_mae']:.3f},"
            f"{summary['delta_eta']:.3f},{summary['delta_pr']:.3f},"
            f"{summary['delta_work']:.3f},{paired}"
        )
    print(f"\nDelta columns are scenario minus {baseline_name}; they are paired operating-point comparisons.")
    print("Pressure/efficiency inputs are separate sources; pressure ratios are interpolated only inside each measured speedline.")
    print("The efficiency digitization is approximate. Ablation effects are model sensitivity, not proof of physical causality.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
