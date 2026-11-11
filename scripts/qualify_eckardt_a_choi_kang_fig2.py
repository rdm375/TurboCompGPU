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


def interpolate_y(points: list[tuple[float, float]], x: float) -> float:
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    raise ValueError("Could not bracket interpolation query")


def crossings(points: list[tuple[float, float]], target: float) -> list[float]:
    result: list[float] = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == target:
            result.append(x0)
        if (y0 - target) * (y1 - target) < 0.0:
            result.append(x0 + (target - y0) * (x1 - x0) / (y1 - y0))
    if points[-1][1] == target:
        result.append(points[-1][0])
    return sorted(set(result))


def uncertainty_intervals(
    points: list[tuple[float, float]], tolerance_pp: float
) -> list[tuple[float, float]]:
    boundaries = [points[0][0], points[-1][0]]
    boundaries.extend(crossings(points, -tolerance_pp))
    boundaries.extend(crossings(points, tolerance_pp))
    boundaries = sorted(set(boundaries))
    intervals: list[tuple[float, float]] = []
    for lower, upper in zip(boundaries, boundaries[1:]):
        middle = 0.5 * (lower + upper)
        if abs(interpolate_y(points, middle)) <= tolerance_pp:
            if intervals and abs(intervals[-1][1] - lower) < 1.0e-10:
                intervals[-1] = (intervals[-1][0], upper)
            else:
                intervals.append((lower, upper))
    for point_x, point_y in points:
        if abs(point_y) <= tolerance_pp and not any(
            lower - 1.0e-10 <= point_x <= upper + 1.0e-10
            for lower, upper in intervals
        ):
            intervals.append((point_x, point_x))
    return sorted(intervals)


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
        if scenario.startswith("flow_mix_sweep_"):
            continue
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

    # Interpolate the dense chi sweep only to expose how the matched-plane
    # residual changes with flow. The source digitization uncertainty is
    # approximately +/- 0.015 efficiency, or +/- 1.5 percentage points.
    # This is a diagnostic consistency check, not coefficient fitting.
    digitization_tolerance_pp = 1.5
    mixing_curves: dict[float, list[tuple[float, float]]] = {}
    for row in ablation_predictions:
        scenario = str(row["scenario"])
        if not scenario.startswith("flow_mix_sweep_chi_"):
            continue
        chi_text = scenario.removeprefix("flow_mix_sweep_chi_").split("_bstar_")[0]
        flow = float(row["mass_flow_kg_s"])
        mixing_curves.setdefault(flow, []).append((
            float(chi_text), float(row["efficiency_error_percentage_points"])))
    for curve in mixing_curves.values():
        curve.sort()

    mixing_diagnostics: list[dict[str, float | str]] = []
    mixing_diagnostic_csv: Path | None = None
    if len(mixing_curves) == len(measurements):
        chi_grid = sorted({chi for curve in mixing_curves.values() for chi, _ in curve})
        shared_scores = [
            (chi, mean([abs(interpolate_y(curve, chi))
                        for curve in mixing_curves.values()]))
            for chi in chi_grid
        ]
        shared_chi, shared_mae_pp = min(shared_scores, key=lambda item: item[1])
        for measurement in measurements:
            flow = measurement["mass_flow_kg_s"]
            curve = mixing_curves[flow]
            zero_crossings = crossings(curve, 0.0)
            tolerance_ranges = uncertainty_intervals(
                curve, digitization_tolerance_pp)
            mixing_diagnostics.append({
                "rpm": RPM,
                "mass_flow_kg_s": flow,
                "measured_efficiency": measurement[
                    "measured_mixed_plane_efficiency"],
                "digitization_tolerance_pp": digitization_tolerance_pp,
                "chi_at_zero_error": zero_crossings[0] if zero_crossings else math.nan,
                "chi_within_uncertainty_min": (
                    tolerance_ranges[0][0] if tolerance_ranges else math.nan),
                "chi_within_uncertainty_max": (
                    tolerance_ranges[-1][1] if tolerance_ranges else math.nan),
                "shared_chi_minimum_mae": shared_chi,
                "efficiency_error_at_shared_chi_pp": interpolate_y(curve, shared_chi),
                "shared_chi_mae_pp": shared_mae_pp,
                "efficiency_error_at_reference_chi_0p15_pp": interpolate_y(curve, 0.15),
            })
        mixing_diagnostic_csv = outdir / "mixing_flow_diagnostic.csv"
        with mixing_diagnostic_csv.open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=list(mixing_diagnostics[0]))
            writer.writeheader()
            writer.writerows(mixing_diagnostics)
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
    if mixing_diagnostics:
        report.extend([
            "",
            "## Flow dependence of the mixing sensitivity",
            "",
            "A dense chi sweep at b*=1.00 was evaluated at the same three measured "
            "flows and station planes. The reported chi values are piecewise-linear "
            "interpolations between model runs. The interval column shows where the "
            "efficiency residual falls within the approximate +/-1.5 percentage-point "
            "digitization uncertainty.",
            "",
            f"The best single chi on the sampled 0.01 grid is {shared_chi:.2f}; its "
            f"three-point efficiency MAE is {shared_mae_pp:.2f} percentage points. "
            "This common-chi score is a sensitivity diagnostic, not an independently "
            "validated calibration.",
            "",
            "| Flow (kg/s) | chi at zero residual | chi interval within digitization uncertainty | "
            "Error at shared chi (pp) | Error at reference chi=0.15 (pp) |",
            "|---:|---:|---:|---:|---:|",
        ])
        for row in mixing_diagnostics:
            chi_zero = row["chi_at_zero_error"]
            chi_lo = row["chi_within_uncertainty_min"]
            chi_hi = row["chi_within_uncertainty_max"]
            chi_text = "none" if math.isnan(float(chi_lo)) else (
                f"{float(chi_lo):.3f}–{float(chi_hi):.3f}")
            zero_text = "not bracketed" if math.isnan(float(chi_zero)) else f"{float(chi_zero):.3f}"
            report.append(
                f"| {row['mass_flow_kg_s']:.3f} | {zero_text} | {chi_text} | "
                f"{row['efficiency_error_at_shared_chi_pp']:+.2f} | "
                f"{row['efficiency_error_at_reference_chi_0p15_pp']:+.2f} |")
        report.extend([
            "",
            "A shared chi is consistent with all three points only if the three "
            "uncertainty intervals overlap. Sparse, approximate digitization cannot "
            "establish a physical flow-dependent mixing law; additional measured "
            "points are needed before changing the closure.",
            "",
            "Pointwise mixing diagnostic: `mixing_flow_diagnostic.csv`.",
        ])
    else:
        report.extend([
            "",
            "## Mixing flow diagnostic",
            "",
            "Not run: the current benchmark source does not emit the dense chi sweep "
            "required by the historical flow-dependence diagnostic. The measured-plane "
            "comparison and the emitted named loss scenarios above are complete.",
        ])

    report_path = outdir / "qualification.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    print(f"\nCSV: {csv_path}")
    print(f"Report: {report_path}")
    print(f"Matched-plane loss ablation: {ablation_csv}")
    print(f"Matched-plane loss summary: {ablation_summary_csv}")
    if mixing_diagnostic_csv is not None:
        print(f"Mixing flow diagnostic: {mixing_diagnostic_csv}")
    else:
        print("Mixing flow diagnostic: not emitted by current benchmark source")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
