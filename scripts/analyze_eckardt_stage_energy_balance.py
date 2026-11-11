#!/usr/bin/env python3
"""Summarize stage energy terms in the Eckardt ablation CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ablation_csv", type=Path)
    parser.add_argument("--rpm", type=float, default=14000.0)
    parser.add_argument(
        "--scenarios",
        nargs="*",
        default=[
            "all_losses_enthalpy_equivalent",
            "published_eckardt_a_set_2042",
        ],
    )
    args = parser.parse_args()

    required = {
        "scenario",
        "rpm",
        "mass_flow_kg_s",
        "stage_total_enthalpy_rise_j_kg",
        "impeller_euler_work_j_kg",
        "external_loss_enthalpy_addition_j_kg",
        "stage_energy_residual_vs_euler_j_kg",
        "stage_energy_residual_vs_euler_disk_recirc_j_kg",
        "stage_energy_residual_vs_euler_external_addition_j_kg",
    }

    with args.ablation_csv.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "CSV lacks energy audit columns: "
                + ", ".join(sorted(missing))
                + "; rebuild and rerun the updated benchmark first"
            )
        rows = [
            row for row in reader
            if abs(float(row["rpm"]) - args.rpm) < 0.5
            and row["scenario"] in args.scenarios
        ]

    if not rows:
        raise ValueError("No rows matched the requested rpm and scenarios")

    print("PR-air stage energy audit; residuals are Δh0 minus named terms")
    print(
        "scenario,flow_kg_s,stage_dh0_kJ_kg,euler_kJ_kg,"
        "external_addition_kJ_kg,resid_vs_euler_kJ_kg,"
        "resid_vs_euler_disk_recirc_kJ_kg,"
        "resid_vs_euler_external_addition_kJ_kg"
    )

    for row in sorted(
        rows,
        key=lambda item: (item["scenario"], float(item["mass_flow_kg_s"])),
    ):
        def kj(column: str) -> float:
            return float(row[column]) / 1000.0

        print(
            f'{row["scenario"]},{float(row["mass_flow_kg_s"]):.3f},'
            f'{kj("stage_total_enthalpy_rise_j_kg"):+.3f},'
            f'{kj("impeller_euler_work_j_kg"):+.3f},'
            f'{kj("external_loss_enthalpy_addition_j_kg"):+.3f},'
            f'{kj("stage_energy_residual_vs_euler_j_kg"):+.3f},'
            f'{kj("stage_energy_residual_vs_euler_disk_recirc_j_kg"):+.3f},'
            f'{kj("stage_energy_residual_vs_euler_external_addition_j_kg"):+.3f}'
        )

    print(
        "These columns report implemented bookkeeping; they do not assume "
        "every equivalent loss term is external shaft work."
    )
    return 0


if __name__ == "__main__":
    main()
