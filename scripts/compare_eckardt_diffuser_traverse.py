#!/usr/bin/env python3
"""Screen radial-discharge Eckardt O traverse data against the Impeller A model."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS = ROOT / "data/validation/eckardt_radial_impeller_o_diffuser_traverse.csv"
DEFAULT_OUTPUT_DIR = ROOT / "benchmarks/results-eckardt-cross-family-screening/radial_o_data_vs_a_model"
NEAR = "near_1.075"
FAR = "far_1.6875"


def read_measurements(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    if not rows:
        raise ValueError(f"No traverse measurement rows in {path}")
    return [{key: float(value) for key, value in row.items()} for row in rows]


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def pressure_metrics(rows: list[dict[str, object]], station: str) -> dict[str, float]:
    subset = [r for r in rows if r["station_id"] == station]
    errors = [float(r["pressure_error_pct"]) for r in subset]
    return {
        "n": len(errors),
        "bias_pct": mean(errors),
        "mape_pct": mean([abs(e) for e in errors]),
        "max_abs_pct": max(abs(e) for e in errors),
    }


def angle_metrics(rows: list[dict[str, object]]) -> dict[str, float]:
    errors = [float(r["angle_error_deg"]) for r in rows if r["angle_error_deg"] != ""]
    return {
        "n": len(errors),
        "bias_deg": mean(errors),
        "mae_deg": mean([abs(e) for e in errors]),
        "max_abs_deg": max(abs(e) for e in errors),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True,
                        help="Built turbocompgpu_eckardt_diffuser_friction_sweep executable")
    parser.add_argument("--measurements", type=Path, default=MEASUREMENTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--flow-basis", choices=("orifice", "traverse"), default="orifice",
                        help="Operating flow supplied to the model (default: orifice measurement)")
    parser.add_argument("--allow-cross-family-screening", action="store_true",
                        help="Allow exploratory Impeller O data vs A model comparison; not A validation")
    args = parser.parse_args()

    if not args.allow_cross_family_screening:
        parser.error("measurements are radial-discharge Impeller O, while this executable uses backswept Impeller A; pass --allow-cross-family-screening only for an exploratory screen")

    exe = args.executable.expanduser().resolve()
    if not exe.is_file():
        parser.error(f"benchmark executable does not exist: {exe}")
    measured = read_measurements(args.measurements)
    flow_column = ("orifice_mass_flow_kg_s" if args.flow_basis == "orifice"
                   else "traverse_mass_flow_kg_s")
    outdir = args.output_dir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # The existing benchmark accepts a pressure-ratio column for map checks;
    # traverse predictions depend on rpm and mass flow only. A dummy 1.0 is
    # supplied here and the benchmark's unrelated map summary is suppressed.
    with tempfile.TemporaryDirectory(prefix="eckardt-traverse-") as temp:
        tempdir = Path(temp)
        operating_points = tempdir / "operating_points.csv"
        with operating_points.open("w", encoding="utf-8") as f:
            f.write("# traverse qualification input; third column is unused by traverse predictions\n")
            for row in measured:
                f.write(f"{row['rpm']:.0f},{row[flow_column]:.8g},1.0\n")
        pressure_output = tempdir / "unused_pressure_summary.csv"
        run = subprocess.run(
            [str(exe), str(operating_points), str(pressure_output)],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            raise RuntimeError(f"Benchmark failed (exit {run.returncode}):\n{run.stdout}")
        prediction_path = tempdir / "impeller_a_diffuser_traverse_predictions.csv"
        with prediction_path.open(newline="", encoding="utf-8") as f:
            predictions = list(csv.DictReader(f))

    if any(r.get("model_geometry") != "eckardt_impeller_a_backswept" for r in predictions):
        raise RuntimeError("Expected backswept Impeller A predictions; refusing a mislabeled comparison")

    by_key = {
        (int(float(r["rpm"])), round(float(r["mass_flow_kg_s"]), 8), r["station_id"]): r
        for r in predictions
    }
    results: list[dict[str, object]] = []
    rotor_static_results: list[dict[str, object]] = []
    for m in measured:
        rpm = int(m["rpm"])
        flow = m[flow_column]
        for station, measured_p, measured_angle in (
            (NEAR, m["near_total_pressure_pa"], m["near_flow_angle_from_meridional_deg"]),
            (FAR, m["far_total_pressure_pa"], None),
        ):
            prediction = by_key[(rpm, round(flow, 8), station)]
            predicted_p = float(prediction["predicted_total_pressure_pa"])
            angle = prediction["predicted_flow_angle_from_meridional_deg"]
            angle_error: float | str = "" if measured_angle is None else float(angle) - measured_angle
            measured_static = None if station == NEAR else m["far_static_pressure_pa"]
            predicted_static = float(prediction["predicted_bulk_static_pressure_pa"])
            static_error: float | str = (
                "" if measured_static is None else
                100.0 * (predicted_static - measured_static) / measured_static)
            measured_p2 = (m["impeller_exit_static_pressure_pa"]
                           if station == NEAR else None)
            predicted_p2 = float(prediction["predicted_impeller_exit_static_pressure_pa"])
            p2_error: float | str = (
                "" if measured_p2 is None else
                100.0 * (predicted_p2 - measured_p2) / measured_p2)
            if measured_p2 is not None:
                rotor_static_results.append({
                    "station_id": "p2",
                    "pressure_error_pct": p2_error,
                })
            results.append({
                "rpm": rpm,
                "orifice_mass_flow_kg_s": m["orifice_mass_flow_kg_s"],
                "traverse_integrated_mass_flow_kg_s": m["traverse_mass_flow_kg_s"],
                "model_mass_flow_kg_s": flow,
                "model_flow_basis": args.flow_basis,
                "station_id": station,
                "r_over_r2": prediction["r_over_r2"],
                "measured_total_pressure_pa": measured_p,
                "predicted_total_pressure_pa": predicted_p,
                "pressure_error_pct": 100.0 * (predicted_p - measured_p) / measured_p,
                "measured_impeller_exit_static_pressure_pa": "" if measured_p2 is None else measured_p2,
                "predicted_impeller_exit_static_pressure_pa": predicted_p2,
                "impeller_exit_static_pressure_error_pct": p2_error,
                "measured_static_pressure_pa": "" if measured_static is None else measured_static,
                "predicted_bulk_static_pressure_pa": predicted_static,
                "static_pressure_error_pct": static_error,
                "measured_flow_angle_from_meridional_deg": "" if measured_angle is None else measured_angle,
                "predicted_flow_angle_from_meridional_deg": angle,
                "angle_error_deg": angle_error,
            })

    csv_path = outdir / f"eckardt_o_data_vs_a_model_{args.flow_basis}_flow.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    near_p = pressure_metrics(results, NEAR)
    far_p = pressure_metrics(results, FAR)
    p2_metrics = pressure_metrics(rotor_static_results, "p2")
    angles = angle_metrics(results)
    measured_losses, predicted_losses = [], []
    far_static_results: list[dict[str, object]] = []
    for m in measured:
        key = (int(m["rpm"]), round(m[flow_column], 8))
        near = by_key[(*key, NEAR)]
        far = by_key[(*key, FAR)]
        measured_losses.append(100.0 * (1.0 - m["far_total_pressure_pa"] / m["near_total_pressure_pa"]))
        predicted_losses.append(100.0 * (1.0 - float(far["predicted_total_pressure_pa"]) /
                                           float(near["predicted_total_pressure_pa"])))
        far_static_results.append({
            "station_id": FAR,
            "pressure_error_pct": 100.0 * (float(far["predicted_bulk_static_pressure_pa"]) -
                                             m["far_static_pressure_pa"]) /
                                  m["far_static_pressure_pa"],
        })
    far_static = pressure_metrics(far_static_results, FAR)

    summary = [
        "# Cross-family screen: Impeller O data vs backswept Impeller A model",
        "",
        "**Not an Impeller A qualification.** The measurements are for radial-discharge Impeller O; the executable uses backswept Impeller A geometry. These residuals cannot be interpreted as Impeller A model errors.",
        "",
        "Compared the model with six published 14,000 rpm mass-averaged points for exploratory cross-family context only.",
        f"Model points are evaluated using the published {args.flow_basis} mass-flow column. Near-station flow angle is compared only at r/r2=1.075; the paper's Table 4 does not report a far-station angle.",
        "",
        "| Quantity | Points | Bias | Mean absolute error | Maximum absolute error |",
        "|---|---:|---:|---:|---:|",
        f"| Impeller-exit static pressure P2 | {p2_metrics['n']} | {p2_metrics['bias_pct']:+.2f}% | {p2_metrics['mape_pct']:.2f}% | {p2_metrics['max_abs_pct']:.2f}% |",
        f"| Near total pressure, r/r2=1.075 | {near_p['n']} | {near_p['bias_pct']:+.2f}% | {near_p['mape_pct']:.2f}% | {near_p['max_abs_pct']:.2f}% |",
        f"| Far total pressure, r/r2=1.6875 | {far_p['n']} | {far_p['bias_pct']:+.2f}% | {far_p['mape_pct']:.2f}% | {far_p['max_abs_pct']:.2f}% |",
        f"| Far static pressure, r/r2=1.6875 | {far_static['n']} | {far_static['bias_pct']:+.2f}% | {far_static['mape_pct']:.2f}% | {far_static['max_abs_pct']:.2f}% |",
        f"| Near flow angle | {angles['n']} | {angles['bias_deg']:+.2f} deg | {angles['mae_deg']:.2f} deg | {angles['max_abs_deg']:.2f} deg |",
        f"| Diffuser total-pressure loss (near to far) | {len(measured)} | measured mean {mean(measured_losses):.2f}%; predicted mean {mean(predicted_losses):.2f}% | mean absolute difference {mean([abs(a-b) for a,b in zip(predicted_losses, measured_losses)]):.2f} percentage points | max absolute difference {max(abs(a-b) for a,b in zip(predicted_losses, measured_losses)):.2f} percentage points |",
        "",
        "The source identifies this single-speed-line dataset as radial-discharge Impeller O and says comparable public traverse data for the backswept impeller was unavailable. Do not use this cross-family screen to tune or qualify Impeller A.",
        "",
        f"Pointwise results: `{csv_path.name}`.",
    ]
    report_path = outdir / f"eckardt_o_data_vs_a_model_{args.flow_basis}_flow.md"
    report_path.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print(f"\nCSV: {csv_path}")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
