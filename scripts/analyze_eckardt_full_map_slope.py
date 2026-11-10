#!/usr/bin/env python3
"""Compare model efficiency levels and flow slopes over Eckardt A digitization."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_efficiency_map(path: Path) -> dict[tuple[int, float], float]:
    rows: dict[tuple[int, float], float] = {}
    with path.open(newline="", encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = next(csv.reader([line]))
            if len(fields) != 3:
                raise ValueError(f"Expected three CSV fields at {path}:{line_number}")
            key = (int(float(fields[0])), round(float(fields[1]), 6))
            rows[key] = float(fields[2])
    if not rows:
        raise ValueError(f"No efficiency points found in {path}")
    return rows


def read_ablation(path: Path, digitized: dict[tuple[int, float], float]):
    grouped = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {
            "scenario", "rpm", "mass_flow_kg_s", "measured_eta_tt",
            "predicted_eta_tt", "eta_error_percentage_points",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            rpm = int(float(row["rpm"]))
            flow = float(row["mass_flow_kg_s"])
            key = (rpm, round(flow, 6))
            if key not in digitized:
                raise ValueError(
                    f"Ablation point {rpm} rpm / {flow:g} kg/s is absent from {path}"
                )
            measured = float(row["measured_eta_tt"])
            if not math.isclose(measured, digitized[key], abs_tol=5.0e-7):
                raise ValueError(
                    f"Stale or mismatched digitization at {rpm} rpm / {flow:g} kg/s: "
                    f"ablation has {measured:g}, map has {digitized[key]:g}"
                )
            predicted = float(row["predicted_eta_tt"])
            error_pp = float(row["eta_error_percentage_points"])
            if not math.isclose(error_pp, 100.0 * (predicted - measured),
                                rel_tol=1.0e-7, abs_tol=1.0e-6):
                raise ValueError(
                    f"Inconsistent efficiency residual at {rpm} rpm / {flow:g} kg/s"
                )
            grouped[row["scenario"]].append({
                "rpm": rpm,
                "flow": flow,
                "measured": measured,
                "predicted": predicted,
                "error_pp": error_pp,
            })
    if not grouped:
        raise ValueError(f"No ablation rows found in {path}")
    return grouped


def slope(rows, field: str) -> float:
    x = [r["flow"] for r in rows]
    y = [r[field] for r in rows]
    xbar = sum(x) / len(x)
    ybar = sum(y) / len(y)
    denominator = sum((value - xbar) ** 2 for value in x)
    if len(rows) < 2 or denominator == 0.0:
        raise ValueError("At least two distinct flows are required for a slope")
    return sum((xv - xbar) * (yv - ybar)
               for xv, yv in zip(x, y)) / denominator


def summarize(rows, scope: str) -> dict[str, object]:
    measured_slope = slope(rows, "measured")
    model_slope = slope(rows, "predicted")
    errors = [r["error_pp"] for r in rows]
    return {
        "scope": scope,
        "rpm": rows[0]["rpm"],
        "points": len(rows),
        "measured_eta_slope_pp_per_kg_s": 100.0 * measured_slope,
        "model_eta_slope_pp_per_kg_s": 100.0 * model_slope,
        "slope_error_model_minus_measured_pp_per_kg_s":
            100.0 * (model_slope - measured_slope),
        "eta_bias_pp": sum(errors) / len(errors),
        "eta_MAE_pp": sum(abs(error) for error in errors) / len(errors),
        "eta_max_abs_pp": max(abs(error) for error in errors),
        "points_within_1pp": sum(abs(error) <= 1.0 for error in errors),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ablation-csv", type=Path, required=True)
    parser.add_argument("--efficiency-map", type=Path,
                        default=Path("data/validation/eckardt_impeller_a_efficiency_digitized.csv"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("benchmarks/results-eckardt-impeller-a/full-map-slope"))
    args = parser.parse_args()

    digitized = read_efficiency_map(args.efficiency_map)
    scenarios = read_ablation(args.ablation_csv, digitized)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pointwise_rows = []
    adjacent_rows = []
    summary_rows = []
    for scenario, rows in sorted(scenarios.items()):
        by_speed = defaultdict(list)
        for row in rows:
            by_speed[row["rpm"]].append(row)
        for rpm, speed_rows in sorted(by_speed.items()):
            speed_rows.sort(key=lambda row: row["flow"])
            for row in speed_rows:
                pointwise_rows.append({
                    "scenario": scenario,
                    "rpm": rpm,
                    "mass_flow_kg_s": row["flow"],
                    "measured_eta_tt": row["measured"],
                    "predicted_eta_tt": row["predicted"],
                    "residual_percentage_points": row["error_pp"],
                })
            for low, high in zip(speed_rows, speed_rows[1:]):
                delta_flow = high["flow"] - low["flow"]
                adjacent_rows.append({
                    "scenario": scenario,
                    "rpm": rpm,
                    "flow_low_kg_s": low["flow"],
                    "flow_high_kg_s": high["flow"],
                    "measured_eta_low": low["measured"],
                    "measured_eta_high": high["measured"],
                    "predicted_eta_low": low["predicted"],
                    "predicted_eta_high": high["predicted"],
                    "measured_secant_slope_pp_per_kg_s":
                        100.0 * (high["measured"] - low["measured"]) / delta_flow,
                    "model_secant_slope_pp_per_kg_s":
                        100.0 * (high["predicted"] - low["predicted"]) / delta_flow,
                    "slope_error_model_minus_measured_pp_per_kg_s":
                        (high["error_pp"] - low["error_pp"]) / delta_flow,
                })
            summary_rows.append({"scenario": scenario,
                                 **summarize(speed_rows, "all_digitized_points")})
            # The 14,000 rpm line contains a low efficiency symbol at 5.69 kg/s.
            # Keep it in the primary result and show its leverage separately.
            if rpm == 14000 and any(abs(row["flow"] - 5.69) < 1.0e-6
                                    for row in speed_rows):
                sensitivity_rows = [row for row in speed_rows
                                    if abs(row["flow"] - 5.69) >= 1.0e-6]
                if len(sensitivity_rows) >= 2:
                    summary_rows.append({
                        "scenario": scenario,
                        **summarize(sensitivity_rows,
                                    "14k_sensitivity_without_5.69kg_s"),
                    })

    pointwise_path = args.output_dir / "pointwise_efficiency_residuals.csv"
    pointwise_fields = ["scenario", "rpm", "mass_flow_kg_s", "measured_eta_tt",
                        "predicted_eta_tt", "residual_percentage_points"]
    with pointwise_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=pointwise_fields)
        writer.writeheader()
        writer.writerows(pointwise_rows)

    adjacent_path = args.output_dir / "adjacent_flow_slopes.csv"
    adjacent_fields = [
        "scenario", "rpm", "flow_low_kg_s", "flow_high_kg_s",
        "measured_eta_low", "measured_eta_high", "predicted_eta_low",
        "predicted_eta_high", "measured_secant_slope_pp_per_kg_s",
        "model_secant_slope_pp_per_kg_s",
        "slope_error_model_minus_measured_pp_per_kg_s",
    ]
    with adjacent_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=adjacent_fields)
        writer.writeheader()
        writer.writerows(adjacent_rows)

    summary_path = args.output_dir / "speedline_slope_summary.csv"
    fields = ["scenario", "scope", "rpm", "points",
              "measured_eta_slope_pp_per_kg_s",
              "model_eta_slope_pp_per_kg_s",
              "slope_error_model_minus_measured_pp_per_kg_s",
              "eta_bias_pp", "eta_MAE_pp", "eta_max_abs_pp",
              "points_within_1pp"]
    with summary_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    report_path = args.output_dir / "speedline_slope_summary.md"
    selected_scenarios = [
        "all_losses_enthalpy_equivalent", "without_skin_friction",
        "mixing_chi_0.05_bstar_1.00", "mixing_chi_0.10_bstar_1.00",
        "mixing_chi_0.15_bstar_1.00", "mixing_chi_0.20_bstar_1.00",
    ]
    report = [
        "# Eckardt A efficiency slope screen across speed lines", "",
        f"Digitized efficiency points: {len(digitized)}. "
        "The plot symbols are approximate (about ±0.01 absolute efficiency).",
        "Ablation rows are checked against the current digitization before analysis.",
        "Slopes are unweighted least-squares fits of efficiency versus mass flow.",
        "The full 14,000 rpm line is the primary result; a second row omits its "
        "5.69 kg/s symbol only to show that point's influence.", "",
        "| Scenario | Scope | rpm | n | Measured slope | Model slope | Slope error | "
        "Efficiency bias | MAE | Max abs | Within ±1 pp |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    selected_rows = [row for row in summary_rows
                     if row["scenario"] in selected_scenarios]
    for row in selected_rows:
        report.append(
            f"| {row['scenario']} | {row['scope']} | {row['rpm']} | {row['points']} | "
            f"{row['measured_eta_slope_pp_per_kg_s']:+.3f} | "
            f"{row['model_eta_slope_pp_per_kg_s']:+.3f} | "
            f"{row['slope_error_model_minus_measured_pp_per_kg_s']:+.3f} | "
            f"{row['eta_bias_pp']:+.2f} pp | {row['eta_MAE_pp']:.2f} pp | "
            f"{row['eta_max_abs_pp']:.2f} pp | {row['points_within_1pp']}/{row['points']} |"
        )
    report.extend([
        "", "## Adjacent-flow slopes for the baseline model", "",
        "These secant slopes compare each neighboring pair of digitized flows. "
        "They show where the residual changes along a speed line and avoid "
        "compressing a local kink into one whole-line regression.", "",
        "| rpm | Flow interval (kg/s) | Measured slope | Model slope | Slope error |",
        "|---:|---:|---:|---:|---:|",
    ])
    baseline_adjacent = [
        row for row in adjacent_rows
        if row["scenario"] == "all_losses_enthalpy_equivalent"
    ]
    for row in baseline_adjacent:
        report.append(
            f"| {row['rpm']} | {row['flow_low_kg_s']:.3f}–{row['flow_high_kg_s']:.3f} | "
            f"{row['measured_secant_slope_pp_per_kg_s']:+.3f} | "
            f"{row['model_secant_slope_pp_per_kg_s']:+.3f} | "
            f"{row['slope_error_model_minus_measured_pp_per_kg_s']:+.3f} |"
        )
    report.extend([
        "", "The skin-friction and chi rows are one-factor sensitivity cases; "
        "they do not establish the physical correctness of either setting.",
        "Agreement on approximate digitized points is a model screen, not "
        "independent validation or a calibrated loss law.", "",
        "Complete summary: `speedline_slope_summary.csv`; point residuals: "
        "`pointwise_efficiency_residuals.csv`; adjacent-pair slopes: "
        "`adjacent_flow_slopes.csv`.",
    ])
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Efficiency points checked: {len(digitized)}")
    print(f"Ablation scenarios: {len(scenarios)}")
    print(f"Speed-line summary: {summary_path}")
    print(f"Pointwise residuals: {pointwise_path}")
    print(f"Adjacent-flow slopes: {adjacent_path}")
    print(f"Report: {report_path}")
    for row in selected_rows:
        print(
            f"{row['scenario']}, {row['scope']}, {row['rpm']} rpm: "
            f"slope measured/model={row['measured_eta_slope_pp_per_kg_s']:+.3f}/"
            f"{row['model_eta_slope_pp_per_kg_s']:+.3f} pp/(kg/s), "
            f"slope error={row['slope_error_model_minus_measured_pp_per_kg_s']:+.3f}, "
            f"eta MAE={row['eta_MAE_pp']:.2f} pp"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
