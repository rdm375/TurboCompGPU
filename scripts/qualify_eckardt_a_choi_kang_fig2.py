#!/usr/bin/env python3
"""Compare Eckardt A predictions with Choi-Kang's three Fig. 2 measurements.

The source defines total pressure at r/r2=1.075, total temperature at
r/r2=1.687, and efficiency from those mixed planes. The benchmark executable
must be built from this checkout so it writes both station predictions.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS = ROOT / "data/validation/eckardt_impeller_a_choi_kang_1999_fig2.csv"
DEFAULT_OUTPUT = ROOT / "benchmarks/results-eckardt-impeller-a/choi-kang-1999-fig2"
RPM = 14000.0
M_OPT = 4.54
P01_PA = 101000.0
T01_K = 288.15
CP_J_KG_K = 1004.5
R_J_KG_K = 287.05
D2_M = 0.4


def read_measurements(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as source:
        for raw in source:
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            row = next(csv.DictReader([raw], fieldnames=[
                "normalized_flow_ratio", "mass_flow_kg_s",
                "measured_total_pressure_ratio_r1_075",
                "measured_cp_dT_over_U2sq_r1_687",
                "measured_mixed_plane_efficiency",
            ]))
            rows.append({key: float(value) for key, value in row.items()})
    if len(rows) != 3:
        raise ValueError(f"Expected three digitized Figure 2 points; got {len(rows)}")
    for row in rows:
        expected = row["normalized_flow_ratio"] * M_OPT
        if abs(row["mass_flow_kg_s"] - expected) > 1.0e-8:
            raise ValueError("Mass-flow values do not match normalized flow and m_opt")
    return rows


def read_station_predictions(path: Path) -> dict[tuple[float, str], dict[str, str]]:
    predictions: dict[tuple[float, str], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            if row["model_geometry"] != "eckardt_impeller_a_backswept":
                raise ValueError(f"Unexpected model geometry: {row['model_geometry']}")
            key = (float(row["mass_flow_kg_s"]), row["station_id"])
            predictions[key] = row
    return predictions


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def read_loss_ablation(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True,
                        help="Built turbocompgpu_eckardt_diffuser_friction_sweep executable")
    parser.add_argument("--measurements", type=Path, default=MEASUREMENTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    executable = args.executable.expanduser().resolve()
    if not executable.is_file():
        parser.error(f"benchmark executable does not exist: {executable}")
    measurements = read_measurements(args.measurements.expanduser().resolve())
    outdir = args.output_dir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="eckardt-a-choi-kang-") as temp:
        tempdir = Path(temp)
        operating_points = tempdir / "operating_points.csv"
        with operating_points.open("w", encoding="utf-8") as target:
            target.write("# The third column is unused by traverse output.\n")
            for row in measurements:
                target.write(f"{RPM:.0f},{row['mass_flow_kg_s']:.8g},1.0\n")
        efficiency_points = tempdir / "efficiency_points.csv"
        with efficiency_points.open("w", encoding="utf-8") as target:
            target.write("# Activate benchmark loss-ablation output at the same flows.\n")
            for row in measurements:
                target.write(
                    f"{RPM:.0f},{row['mass_flow_kg_s']:.8g},"
                    f"{row['measured_mixed_plane_efficiency']:.8g}\n")
        summary = tempdir / "unused_map_summary.csv"
        efficiency_output = tempdir / "unused_efficiency_output.csv"
        run = subprocess.run(
            [str(executable), str(operating_points), str(summary),
             str(efficiency_points), str(efficiency_output)],
            check=False, text=True, capture_output=True,
        )
        if run.returncode:
            raise RuntimeError(
                f"Benchmark failed ({run.returncode}):\n{run.stdout}\n{run.stderr}")
        station_file = tempdir / "impeller_a_diffuser_traverse_predictions.csv"
        if not station_file.is_file():
            raise RuntimeError(
                "Benchmark did not write impeller_a_diffuser_traverse_predictions.csv; "
                "rebuild it from the current checkout.")
        predictions = read_station_predictions(station_file)
        ablation_path = tempdir / "loss_term_ablation.csv"
        if not ablation_path.is_file():
            raise RuntimeError(
                "Benchmark did not write loss_term_ablation.csv; rebuild it "
                "from the current checkout.")
        ablation_rows = read_loss_ablation(ablation_path)

    gamma = CP_J_KG_K / (CP_J_KG_K - R_J_KG_K)
    u2 = math.pi * D2_M * RPM / 60.0
    output_rows: list[dict[str, float]] = []
    for measured in measurements:
        flow = measured["mass_flow_kg_s"]
        near = predictions[(flow, "near_1.075")]
        far = predictions[(flow, "far_1.6875")]
        near_pressure_ratio = float(near["predicted_total_pressure_pa"]) / P01_PA
        far_total_temperature = float(far["predicted_total_temperature_k"])
        temperature_rise_coefficient = (
            CP_J_KG_K * (far_total_temperature - T01_K) / (u2 * u2))
        isentropic_temperature_rise = T01_K * (
            near_pressure_ratio ** ((gamma - 1.0) / gamma) - 1.0)
        actual_temperature_rise = far_total_temperature - T01_K
        efficiency = isentropic_temperature_rise / actual_temperature_rise
        row = dict(measured)
        row.update({
            "predicted_total_pressure_ratio_r1_075": near_pressure_ratio,
            "pressure_ratio_error_percent": 100.0 * (
                near_pressure_ratio / measured["measured_total_pressure_ratio_r1_075"] - 1.0),
            "predicted_cp_dT_over_U2sq_r1_687": temperature_rise_coefficient,
            "temperature_rise_coefficient_error": temperature_rise_coefficient -
                measured["measured_cp_dT_over_U2sq_r1_687"],
            "predicted_mixed_plane_efficiency": efficiency,
            "efficiency_error_percentage_points": 100.0 * (
                efficiency - measured["measured_mixed_plane_efficiency"]),
            "predicted_far_total_temperature_k": far_total_temperature,
        })
        output_rows.append(row)

    # Recompute efficiency from the exact two planes used by Choi and Kang for
    # every ablation scenario. The benchmark's default eta_tt uses its stage
    # inlet/outlet planes, so it is not the comparable response for this screen.
    ablation_predictions: list[dict[str, float | str]] = []
    for raw in ablation_rows:
        flow = float(raw["mass_flow_kg_s"])
        near_pressure_ratio = float(raw["predicted_near_total_pressure_pa"]) / P01_PA
        far_temperature = float(raw["predicted_far_total_temperature_k"])
        efficiency = (
            T01_K * (near_pressure_ratio ** ((gamma - 1.0) / gamma) - 1.0)
            / (far_temperature - T01_K)
        )
        measured = next(
            item for item in measurements
            if abs(item["mass_flow_kg_s"] - flow) < 1.0e-8
        )
        ablation_predictions.append({
            "scenario": raw["scenario"],
            "rpm": float(raw["rpm"]),
            "mass_flow_kg_s": flow,
            "measured_mixed_plane_efficiency": measured[
                "measured_mixed_plane_efficiency"],
            "predicted_near_total_pressure_ratio": near_pressure_ratio,
            "predicted_far_total_temperature_k": far_temperature,
            "predicted_temperature_rise_k": far_temperature - T01_K,
            "predicted_mixed_plane_efficiency": efficiency,
            "efficiency_error_percentage_points": 100.0 * (
                efficiency - measured["measured_mixed_plane_efficiency"]),
        })

    ablation_summaries: list[dict[str, float | str]] = []
    for scenario in dict.fromkeys(str(row["scenario"]) for row in ablation_predictions):
        rows = [row for row in ablation_predictions if row["scenario"] == scenario]
        errors = [float(row["efficiency_error_percentage_points"]) for row in rows]
        ablation_summaries.append({
            "scenario": scenario,
            "points": len(rows),
            "efficiency_bias_pp": mean(errors),
            "efficiency_mae_pp": mean([abs(value) for value in errors]),
            "efficiency_max_abs_pp": max(abs(value) for value in errors),
            "mean_temperature_rise_k": mean([
                float(row["predicted_temperature_rise_k"]) for row in rows]),
            "mean_near_pressure_ratio": mean([
                float(row["predicted_near_total_pressure_ratio"]) for row in rows]),
        })

    ablation_csv = outdir / "station_matched_loss_ablation.csv"
    with ablation_csv.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=list(ablation_predictions[0]))
        writer.writeheader()
        writer.writerows(ablation_predictions)
    ablation_summary_csv = outdir / "station_matched_loss_ablation_summary.csv"
    with ablation_summary_csv.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=list(ablation_summaries[0]))
        writer.writeheader()
        writer.writerows(ablation_summaries)

    csv_path = outdir / "qualification.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    pr_errors = [row["pressure_ratio_error_percent"] for row in output_rows]
    temp_errors = [row["temperature_rise_coefficient_error"] for row in output_rows]
    eta_errors = [row["efficiency_error_percentage_points"] for row in output_rows]
    report = [
        "# Eckardt A / Choi–Kang Figure 2 comparison",
        "",
        "Three approximate experimental symbols digitized from Choi and Kang (1999), Figure 2. "
        "The test is at 14,000 rpm; mass flow is normalized by the reported design flow "
        "4.54 kg/s. The source uses total pressure at r/r2=1.075 and total temperature "
        "at r/r2=1.687 to form its efficiency. Predictions use those same two model stations.",
        "",
        "Digitization uncertainty is estimated at about ±0.015 in pressure ratio, "
        "temperature-rise coefficient, and efficiency. These three points are a "
        "component/performance screen, not broad map qualification.",
        "",
        "| Flow / design flow | Flow (kg/s) | Measured PR | Predicted PR | PR error | "
        "Measured temperature coefficient | Predicted coefficient | "
        "Measured efficiency | Predicted efficiency | Efficiency error |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in output_rows:
        report.append(
            f"| {row['normalized_flow_ratio']:.2f} | {row['mass_flow_kg_s']:.3f} | "
            f"{row['measured_total_pressure_ratio_r1_075']:.3f} | "
            f"{row['predicted_total_pressure_ratio_r1_075']:.4f} | "
            f"{row['pressure_ratio_error_percent']:+.2f}% | "
            f"{row['measured_cp_dT_over_U2sq_r1_687']:.3f} | "
            f"{row['predicted_cp_dT_over_U2sq_r1_687']:.4f} | "
            f"{row['measured_mixed_plane_efficiency']:.3f} | "
            f"{row['predicted_mixed_plane_efficiency']:.4f} | "
            f"{row['efficiency_error_percentage_points']:+.2f} pp |")
    report.extend([
        "",
        f"Across three points: pressure-ratio MAPE {mean([abs(v) for v in pr_errors]):.2f}%; "
        f"mean absolute temperature-coefficient error {mean([abs(v) for v in temp_errors]):.4f}; "
        f"efficiency MAE {mean([abs(v) for v in eta_errors]):.2f} percentage points.",
        "",
        "The source averages experimental quantities from different radial planes; the model "
        "comparison follows that convention. The data are approximate plot digitizations, "
        "and the report keeps their uncertainty visible rather than treating them as exact points.",
        "",
        "Source: Y. S. Choi and S. H. Kang, ‘Application of Through-Flow Calculation to "
        "Design and Performance Prediction of Centrifugal Compressor,’ International Journal "
        "of Rotating Machinery, 5(1), 17–33 (1999), Fig. 2, DOI "
        "[10.1155/S1023621X99000020](https://doi.org/10.1155/S1023621X99000020).",
        "",
        "Pointwise results: `qualification.csv`.",
        "",
        "## Loss-term sensitivity at matched measurement planes",
        "",
        "Each scenario recomputes pressure at r/r2=1.075 and total temperature at "
        "r/r2=1.6875, then forms mixed-plane efficiency using those same stations. "
        "Rows are sensitivity experiments, not evidence that removing a loss term "
        "is physically correct.",
        "",
        "| Scenario | Points | Efficiency bias (pp) | Efficiency MAE (pp) | "
        "Max abs (pp) | Mean ΔT (K) | Mean near PR |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in ablation_summaries:
        report.append(
            f"| {row['scenario']} | {row['points']} | "
            f"{row['efficiency_bias_pp']:+.2f} | {row['efficiency_mae_pp']:.2f} | "
            f"{row['efficiency_max_abs_pp']:.2f} | "
            f"{row['mean_temperature_rise_k']:.2f} | "
            f"{row['mean_near_pressure_ratio']:.4f} |")
    report.extend([
        "",
        "Pointwise ablation: `station_matched_loss_ablation.csv`; summary: "
        "`station_matched_loss_ablation_summary.csv`.",
    ])
    report_path = outdir / "qualification.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    print(f"\nCSV: {csv_path}")
    print(f"Report: {report_path}")
    print(f"Matched-plane loss ablation: {ablation_csv}")
    print(f"Matched-plane loss summary: {ablation_summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
