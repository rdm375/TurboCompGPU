#!/usr/bin/env python3
"""Analyze dense mixing-chi sensitivity against Eckardt A efficiency symbols.

The efficiency points are approximate Figure 6(d) digitizations. They are not
paired pointwise with the separate pressure-ratio transcription, so this is a
same-rotor sensitivity screen rather than an independent full-stage validation.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "data/validation/eckardt_impeller_a_efficiency_digitized.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(row for row in source if row.strip() and not row.startswith("#")))


def read_efficiency_data(path: Path) -> list[dict[str, str]]:
    points: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as source:
        for line in source:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = next(csv.reader([line]))
            if len(fields) != 3:
                raise ValueError(f"Expected 3 efficiency columns, got: {line.rstrip()}")
            points.append({
                "rpm": fields[0],
                "mass_flow_kg_s": fields[1],
                "measured_total_to_total_efficiency": fields[2],
            })
    return points


def key(row: dict[str, str]) -> tuple[int, float]:
    return int(round(float(row["rpm"]))), round(float(row["mass_flow_kg_s"]), 6)


def interpolate(points: list[tuple[float, float]], x: float) -> float:
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    raise ValueError("Could not bracket interpolation query")


def crossings(points: list[tuple[float, float]], target: float) -> list[float]:
    roots: list[float] = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == target:
            roots.append(x0)
        if (y0 - target) * (y1 - target) < 0.0:
            roots.append(x0 + (target - y0) * (x1 - x0) / (y1 - y0))
    if points[-1][1] == target:
        roots.append(points[-1][0])
    return sorted(set(roots))


def tolerance_ranges(
    points: list[tuple[float, float]], tolerance_pp: float
) -> list[tuple[float, float]]:
    bounds = sorted(set([
        points[0][0], points[-1][0],
        *crossings(points, -tolerance_pp),
        *crossings(points, tolerance_pp),
    ]))
    ranges: list[tuple[float, float]] = []
    for low, high in zip(bounds, bounds[1:]):
        if abs(interpolate(points, 0.5 * (low + high))) <= tolerance_pp:
            if ranges and abs(ranges[-1][1] - low) < 1e-10:
                ranges[-1] = (ranges[-1][0], high)
            else:
                ranges.append((low, high))
    for x, y in points:
        if abs(y) <= tolerance_pp and not any(a - 1e-10 <= x <= b + 1e-10 for a, b in ranges):
            ranges.append((x, x))
    return sorted(ranges)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ablation_csv", type=Path,
                        help="loss_term_ablation.csv generated with the dense flow_mix_sweep scenarios")
    parser.add_argument("--efficiency-data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--uncertainty-pp", type=float, default=1.0,
                        help="Approximate absolute efficiency digitization uncertainty in percentage points")
    args = parser.parse_args()
    if not (math.isfinite(args.uncertainty_pp) and args.uncertainty_pp > 0.0):
        parser.error("--uncertainty-pp must be positive and finite")

    measured = read_efficiency_data(args.efficiency_data.expanduser().resolve())
    measured_by_key = {key(row): float(row["measured_total_to_total_efficiency"])
                       for row in measured}
    raw = read_csv(args.ablation_csv.expanduser().resolve())
    curves: dict[tuple[int, float], list[tuple[float, float]]] = {}
    for row in raw:
        name = row["scenario"]
        prefix = "flow_mix_sweep_chi_"
        if not name.startswith(prefix):
            continue
        chi = float(name.removeprefix(prefix).split("_bstar_")[0])
        point_key = key(row)
        if point_key not in measured_by_key:
            continue
        curves.setdefault(point_key, []).append(
            (chi, float(row["eta_error_percentage_points"])))
    for curve in curves.values():
        curve.sort()

    missing = sorted(set(measured_by_key) - set(curves))
    if missing:
        raise SystemExit(
            f"Missing dense chi sweep for {len(missing)} measured points; "
            "rebuild/run the benchmark with the current source")

    chi_grid = sorted({chi for curve in curves.values() for chi, _ in curve})
    scores = [(chi, sum(abs(interpolate(curves[k], chi)) for k in curves) / len(curves))
              for chi in chi_grid]
    shared_chi, shared_mae = min(scores, key=lambda item: item[1])

    output_rows: list[dict[str, object]] = []
    for point_key, measured_eta in sorted(measured_by_key.items()):
        curve = curves[point_key]
        roots = crossings(curve, 0.0)
        ranges = tolerance_ranges(curve, args.uncertainty_pp)
        output_rows.append({
            "rpm": point_key[0],
            "mass_flow_kg_s": point_key[1],
            "measured_eta_tt": measured_eta,
            "chi_at_zero_error": roots[0] if roots else "",
            "chi_within_uncertainty_min": ranges[0][0] if ranges else "",
            "chi_within_uncertainty_max": ranges[-1][1] if ranges else "",
            "error_at_shared_chi_pp": interpolate(curve, shared_chi),
            "error_at_reference_chi_0p15_pp": interpolate(curve, 0.15),
            "uncertainty_pp": args.uncertainty_pp,
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "mixing_flow_diagnostic.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    summaries = []
    for rpm in sorted({row["rpm"] for row in output_rows}):
        rows = [row for row in output_rows if row["rpm"] == rpm]
        errors = [float(row["error_at_shared_chi_pp"]) for row in rows]
        coverage = sum(abs(error) <= args.uncertainty_pp for error in errors)
        summaries.append((rpm, len(rows), sum(errors) / len(errors),
                          sum(abs(error) for error in errors) / len(errors),
                          max(abs(error) for error in errors), coverage))

    pointwise_ranges = [
        (float(row["chi_within_uncertainty_min"]),
         float(row["chi_within_uncertainty_max"]))
        for row in output_rows
        if row["chi_within_uncertainty_min"] != ""
        and row["chi_within_uncertainty_max"] != ""
    ]
    shared_range_text = "no all-point overlap"
    if len(pointwise_ranges) == len(output_rows):
        common_low = max(low for low, _ in pointwise_ranges)
        common_high = min(high for _, high in pointwise_ranges)
        if common_low <= common_high:
            shared_range_text = f"{common_low:.3f}–{common_high:.3f}"

    report = [
        "# Eckardt A mixing-chi flow-dependence screen",
        "",
        f"Efficiency points: {len(measured)} across {len(summaries)} speed lines. "
        "The efficiency symbols are approximate Figure 6(d) digitizations, with "
        "about ±0.01 absolute efficiency uncertainty. They are not pointwise paired "
        "with the separate pressure-ratio transcription.",
        "",
        "The dense mixing sweep uses b*=1.00 and changes chi while retaining the "
        "other benchmark loss settings. Per-point chi intervals are linearly "
        "interpolated between model evaluations. This is a sensitivity screen, "
        "not a fitted or validated flow-dependent closure.",
        "",
        f"Best shared chi on the sampled 0.01 grid: {shared_chi:.2f}; "
        f"{len(output_rows)}-point efficiency MAE: {shared_mae:.2f} percentage points; "
        f"points within ±{args.uncertainty_pp:.2f} pp: "
        f"{sum(abs(float(r['error_at_shared_chi_pp'])) <= args.uncertainty_pp for r in output_rows)}/{len(output_rows)}.",
        f"Intersection of per-point chi intervals within digitization uncertainty: {shared_range_text}.",
        "",
        "| Speed (rpm) | Points | Shared-chi bias (pp) | Shared-chi MAE (pp) | "
        "Max absolute error (pp) | Within uncertainty |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for rpm, count, bias, mae, maximum, within in summaries:
        report.append(f"| {rpm} | {count} | {bias:+.2f} | {mae:.2f} | {maximum:.2f} | {within}/{count} |")
    report.extend([
        "",
        "A single chi is consistent with all digitized points only when their "
        "pointwise uncertainty intervals share a nonempty intersection. Even then, "
        "agreement would be calibration sensitivity, not physical validation, "
        "because the efficiency digitization is approximate and is unpaired with "
        "the pressure-ratio measurements.",
        "",
        "Pointwise output: `mixing_flow_diagnostic.csv`.",
    ])
    report_path = args.output_dir / "mixing_flow_diagnostic.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    print(f"\nCSV: {csv_path}")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
