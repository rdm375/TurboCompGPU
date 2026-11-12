#!/usr/bin/env python3
"""Summarize digitized experimental trends from Rossbach et al. Figure 9."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data/validation/rossbach_etal_2015_diffusion_pressure_recovery_digitized.csv"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    series: dict[str, list[tuple[float, float]]] = {}
    with args.csv.open(newline="") as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            series.setdefault(row["component"], []).append(
                (float(row["phi_stage"]), float(row["measured_cp"]))
            )
    required = {"plane_2_to_4", "plane_2_to_out", "plane_4_to_out"}
    if not required <= series.keys():
        raise SystemExit(f"Missing required series: {sorted(required - series.keys())}")
    for name, points in sorted(series.items()):
        points.sort()
        print(f"{name}: n={len(points)}, Cp range=[{min(y for _, y in points):+.3f}, {max(y for _, y in points):+.3f}]")
    diffuser = series["plane_2_to_4"]
    stage = series["plane_2_to_out"]
    ret = series["plane_4_to_out"]
    diffuser_span = max(y for _, y in diffuser) - min(y for _, y in diffuser)
    stage_delta = stage[-1][1] - stage[0][1]
    high_flow_negative_return = ret[-1][1] < 0.0
    print(f"diffuser Cp span across speedline: {diffuser_span:.3f}")
    print(f"stage Cp change low-to-high flow: {stage_delta:+.3f}")
    print(f"high-flow return-channel Cp negative: {high_flow_negative_return}")
    if not (diffuser_span < abs(stage_delta) and stage_delta < 0 and high_flow_negative_return):
        raise SystemExit("Digitized trends do not match the paper's stated qualitative results; inspect the read-off")
    print("Trend checks agree with the paper's qualitative Figure 9 discussion.")
    print("These are approximate figure-read trends, not a TurboCompGPU validation score.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
