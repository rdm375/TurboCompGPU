#!/usr/bin/env python3
"""Score 14,000 rpm Eckardt A scenarios on efficiency, PR, and work residuals.

The 5.69 kg/s point is reported separately. No weighted scalar score is used;
the ordinary-point candidates are compared by Pareto dominance.
"""
import argparse
import csv
import math
from pathlib import Path

RPM = 14000.0
ORDINARY_FLOWS = (3.99, 4.56, 5.11)
HIGH_FLOW = 5.69
METRICS = ("eta_mae_pp", "pr_mape_percent", "work_mae_kj_kg")
REQUIRED = {
    "model", "rpm", "mass_flow_kg_s", "eta_error_percentage_points",
    "pressure_ratio_pairing", "predicted_pressure_ratio",
    "measured_pressure_ratio_interpolated",
    "specific_work_residual_predicted_minus_implied_J_kg",
    "work_residual_component_from_pressure_ratio_J_kg",
    "work_residual_component_from_efficiency_J_kg",
}


def read_rows(path):
    grouped = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(
            line for line in stream if not line.lstrip().startswith("#")
        )
        missing = REQUIRED.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if not math.isclose(float(row["rpm"]), RPM, abs_tol=1e-8):
                continue
            if row["pressure_ratio_pairing"] != "linear_interpolation":
                continue
            flow = float(row["mass_flow_kg_s"])
            grouped.setdefault(row["model"], {})[flow] = {
                "eta_error_pp": float(row["eta_error_percentage_points"]),
                "predicted_pr": float(row["predicted_pressure_ratio"]),
                "measured_pr": float(row["measured_pressure_ratio_interpolated"]),
                "work_j_kg": float(
                    row["specific_work_residual_predicted_minus_implied_J_kg"]
                ),
                "pr_work_j_kg": float(
                    row["work_residual_component_from_pressure_ratio_J_kg"]
                ),
                "eta_work_j_kg": float(
                    row["work_residual_component_from_efficiency_J_kg"]
                ),
            }
    return grouped


def summarize(name, rows):
    ordinary = [rows[flow] for flow in ORDINARY_FLOWS]
    high = rows[HIGH_FLOW]
    n = len(ordinary)
    return {
        "scenario": name,
        "ordinary_points": n,
        "eta_bias_pp": sum(r["eta_error_pp"] for r in ordinary) / n,
        "eta_mae_pp": sum(abs(r["eta_error_pp"]) for r in ordinary) / n,
        "pr_mape_percent": sum(
            abs(100 * (r["predicted_pr"] / r["measured_pr"] - 1))
            for r in ordinary
        ) / n,
        "work_bias_kj_kg": sum(r["work_j_kg"] for r in ordinary) / n / 1000,
        "work_mae_kj_kg": sum(abs(r["work_j_kg"]) for r in ordinary) / n / 1000,
        "ordinary_pr_work_component_kj_kg":
            sum(r["pr_work_j_kg"] for r in ordinary) / n / 1000,
        "ordinary_eta_work_component_kj_kg":
            sum(r["eta_work_j_kg"] for r in ordinary) / n / 1000,
        "high_flow_eta_error_pp": high["eta_error_pp"],
        "high_flow_pr_error_percent":
            100 * (high["predicted_pr"] / high["measured_pr"] - 1),
        "high_flow_work_residual_kj_kg": high["work_j_kg"] / 1000,
        "high_flow_pr_work_component_kj_kg": high["pr_work_j_kg"] / 1000,
        "high_flow_eta_work_component_kj_kg": high["eta_work_j_kg"] / 1000,
    }


def pareto_front(rows):
    frontier = set()
    for candidate in rows:
        dominated = False
        for other in rows:
            if other is candidate:
                continue
            no_worse = all(other[key] <= candidate[key] for key in METRICS)
            strictly_better = any(other[key] < candidate[key] for key in METRICS)
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            frontier.add(candidate["scenario"])
    return frontier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("energy_audit", type=Path)
    parser.add_argument("--output", type=Path,
                        default=Path("eckardt_joint_residual_scores.csv"))
    args = parser.parse_args()

    grouped = read_rows(args.energy_audit)
    flows = (*ORDINARY_FLOWS, HIGH_FLOW)
    complete = {
        name: rows for name, rows in grouped.items()
        if all(flow in rows for flow in flows)
    }
    if not complete:
        raise ValueError(
            "No 14,000 rpm scenarios contain all four requested flow points"
        )

    summaries = [summarize(name, rows) for name, rows in complete.items()]
    frontier = pareto_front(summaries)
    for row in summaries:
        row["ordinary_pareto"] = (
            "yes" if row["scenario"] in frontier else "no"
        )
    summaries.sort(
        key=lambda row: (
            row["ordinary_pareto"] != "yes",
            row["eta_mae_pp"],
        )
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    print("Ordinary points: 3.99, 4.56, 5.11 kg/s; no weighted score")
    print(
        "scenario,pareto,eta_MAE_pp,PR_MAPE_pct,work_MAE_kJ_kg,"
        "work_bias_kJ_kg,PR_work_component,eta_work_component,"
        "5.69_eta_error_pp,5.69_PR_error_pct,5.69_work_residual_kJ_kg"
    )
    for row in summaries:
        print(
            f"{row['scenario']},{row['ordinary_pareto']},"
            f"{row['eta_mae_pp']:.3f},{row['pr_mape_percent']:.3f},"
            f"{row['work_mae_kj_kg']:.3f},{row['work_bias_kj_kg']:+.3f},"
            f"{row['ordinary_pr_work_component_kj_kg']:+.3f},"
            f"{row['ordinary_eta_work_component_kj_kg']:+.3f},"
            f"{row['high_flow_eta_error_pp']:+.2f},"
            f"{row['high_flow_pr_error_percent']:+.2f},"
            f"{row['high_flow_work_residual_kj_kg']:+.2f}"
        )
    print(f"CSV: {args.output}")
    print("Pareto means no other tested scenario is as good or better on all "
          "three ordinary-point metrics and strictly better on at least one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
