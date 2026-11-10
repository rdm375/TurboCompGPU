#!/usr/bin/env python3
"""Audit the flow dependence of Eckardt A efficiency and component losses.

Reads the detailed loss_term_ablation.csv written by
benchmarks/eckardt_diffuser_friction_sweep.cpp. By default it analyzes the
14 krpm corrected-digitization speed line and omits the disputed 5.69 kg/s
efficiency symbol. The output is a diagnostic trend comparison, not a fit or
proof of physical causality.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SCENARIO = "all_losses_enthalpy_equivalent"
LOSS_COLUMNS = (
    ("skin friction", "impeller_skin_friction_loss_j_kg"),
    ("blade loading", "impeller_blade_loading_loss_j_kg"),
    ("incidence", "impeller_incidence_loss_j_kg"),
    ("clearance", "impeller_clearance_loss_j_kg"),
    ("mixing", "impeller_mixing_loss_j_kg"),
    ("disk friction", "impeller_disk_friction_loss_j_kg"),
    ("diffuser loss", "diffuser_loss_percent"),
)


def slope(xs: list[float], ys: list[float]) -> float:
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    denominator = sum((x - xbar) ** 2 for x in xs)
    if denominator == 0.0:
        raise ValueError("Need at least two distinct mass-flow values")
    return sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denominator


def read_rows(path: Path, rpm: float, excluded_flows: list[float], tolerance: float):
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {
            "scenario", "rpm", "mass_flow_kg_s", "measured_eta_tt",
            "predicted_eta_tt", "eta_error_percentage_points",
            "predicted_total_pressure_ratio", "impeller_specific_work_j_kg",
            *(column for _, column in LOSS_COLUMNS),
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("CSV is missing columns: " + ", ".join(sorted(missing)))
        rows = []
        for row in reader:
            if row["scenario"] != SCENARIO:
                continue
            if abs(float(row["rpm"]) - rpm) > tolerance:
                continue
            flow = float(row["mass_flow_kg_s"])
            if any(abs(flow - excluded) <= tolerance for excluded in excluded_flows):
                continue
            rows.append(row)
    rows.sort(key=lambda row: float(row["mass_flow_kg_s"]))
    if len(rows) < 2:
        raise ValueError(f"Need at least two {rpm:g} rpm points after exclusions")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "ablation_csv", nargs="?", type=Path,
        default=Path("benchmarks/results-eckardt-impeller-a/corrected-digitization/loss_term_ablation.csv"),
        help="Detailed loss_term_ablation.csv from the corrected digitization run",
    )
    parser.add_argument("--rpm", type=float, default=14000.0)
    parser.add_argument(
        "--exclude-flow", type=float, action="append", default=[5.69],
        help="Exclude a mass flow in kg/s (repeatable; default excludes 5.69)",
    )
    parser.add_argument("--tolerance", type=float, default=0.02,
                        help="RPM/flow matching tolerance; default 0.02")
    parser.add_argument("--output-csv", type=Path,
                        help="Optional path for the selected pointwise rows")
    args = parser.parse_args()

    source_path = args.ablation_csv.expanduser().resolve()
    rows = read_rows(source_path, args.rpm, args.exclude_flow, args.tolerance)
    flows = [float(row["mass_flow_kg_s"]) for row in rows]
    measured = [float(row["measured_eta_tt"]) for row in rows]
    predicted = [float(row["predicted_eta_tt"]) for row in rows]
    errors_pp = [100.0 * (p - m) for p, m in zip(predicted, measured)]

    print(f"Source: {source_path}")
    print(f"Scenario: {SCENARIO}; speed: {args.rpm:g} rpm; points: {len(rows)}")
    print("Excluded flows (kg/s): " + ", ".join(f"{x:g}" for x in args.exclude_flow))
    print(f"Measured efficiency slope: {slope(flows, measured):+.5f} per (kg/s)")
    print(f"Model efficiency slope:    {slope(flows, predicted):+.5f} per (kg/s)")
    print(f"Residual slope:            {slope(flows, errors_pp):+.3f} pp per (kg/s)")
    print("\nPointwise values:")
    print("flow_kg_s, eta_measured, eta_model, error_pp, PR_model, work_kJ_kg")
    for row, error in zip(rows, errors_pp):
        print(
            f"{float(row['mass_flow_kg_s']):.3f}, "
            f"{float(row['measured_eta_tt']):.4f}, "
            f"{float(row['predicted_eta_tt']):.4f}, "
            f"{error:+.2f}, "
            f"{float(row['predicted_total_pressure_ratio']):.4f}, "
            f"{float(row['impeller_specific_work_j_kg']) / 1000.0:.3f}"
        )

    print("\nLoss trends (least-squares slope over the selected points):")
    print("term, first_kJ_kg_or_pct, last_kJ_kg_or_pct, slope_per_kg_s")
    for label, column in LOSS_COLUMNS:
        values = [float(row[column]) for row in rows]
        scale = 1000.0 if column.endswith("_j_kg") else 1.0
        y = [value / scale for value in values]
        units = "kJ/kg" if scale == 1000.0 else "%"
        print(
            f"{label} ({units}), {y[0]:.5f}, {y[-1]:.5f}, "
            f"{slope(flows, y):+.5f}"
        )

    if args.output_csv:
        output = args.output_csv.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nSelected rows: {output}")

    print("\nInterpretation limit: these are model loss trends against an approximate, separately digitized efficiency curve. A component's correlation with the residual does not establish causality.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
