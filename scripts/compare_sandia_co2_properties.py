#!/usr/bin/env python3
"""Compare the Sandia PR implementation with CoolProp's CO2 HEOS reference.

This is a diagnostic comparison with no pass/fail accuracy threshold. It runs
the repository's C++ Peng-Robinson implementation to export states, then
compares those states with CoolProp/HEOS (the Span-Wagner CO2 formulation).
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

from CoolProp.CoolProp import PropsSI
import CoolProp


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sandia-co2-properties-") as temp:
        property_csv = Path(temp) / "pr_properties.csv"
        screen_csv = Path(temp) / "sandia_screen.csv"
        subprocess.run(
            [str(args.executable), str(args.map), str(screen_csv), str(property_csv)],
            check=True,
        )

        with property_csv.open(newline="", encoding="utf-8") as source:
            pr_rows = list(csv.DictReader(source))
    if not pr_rows:
        raise RuntimeError("PR executable produced no property rows")

    anchor = next((row for row in pr_rows if row["label"] == "datum_anchor"), None)
    if anchor is None or anchor["status"] != "OK":
        raise RuntimeError("PR executable did not produce the enthalpy/entropy datum anchor")
    T_anchor = float(anchor["T_K"])
    p_anchor = float(anchor["p_Pa"])
    h_anchor_ref = PropsSI("H", "T", T_anchor, "P", p_anchor, "HEOS::CO2")
    s_anchor_ref = PropsSI("S", "T", T_anchor, "P", p_anchor, "HEOS::CO2")
    h_anchor_pr = float(anchor["h_J_kg"])
    s_anchor_pr = float(anchor["s_J_kgK"])

    output_rows: list[dict[str, str | float]] = []
    by_group: dict[str, list[dict[str, str | float]]] = {
        "Sandia operating points": [],
        "near-critical diagnostic grid": [],
        "datum anchor": [],
    }
    for row in pr_rows:
        if row["status"] != "OK":
            continue
        T, p = float(row["T_K"]), float(row["p_Pa"])
        rho_ref = PropsSI("D", "T", T, "P", p, "HEOS::CO2")
        cp_ref = PropsSI("C", "T", T, "P", p, "HEOS::CO2")
        h_ref = PropsSI("H", "T", T, "P", p, "HEOS::CO2")
        s_ref = PropsSI("S", "T", T, "P", p, "HEOS::CO2")
        rho_pr, cp_pr = float(row["rho_kg_m3"]), float(row["cp_J_kgK"])
        delta_h_pr = float(row["h_J_kg"]) - h_anchor_pr
        delta_h_ref = h_ref - h_anchor_ref
        delta_s_pr = float(row["s_J_kgK"]) - s_anchor_pr
        delta_s_ref = s_ref - s_anchor_ref
        rho_error = 100.0 * (rho_pr / rho_ref - 1.0)
        cp_error = 100.0 * (cp_pr / cp_ref - 1.0)
        dh_error = delta_h_pr - delta_h_ref
        ds_error = delta_s_pr - delta_s_ref
        if row["label"].startswith("sandia_rpm_"):
            group = "Sandia operating points"
        elif row["label"].startswith("grid_"):
            group = "near-critical diagnostic grid"
        else:
            group = "datum anchor"
        comparison = {
            "sample_group": group,
            "label": row["label"], "T_K": T, "p_Pa": p,
            "rho_PR_kg_m3": rho_pr, "rho_HEOS_kg_m3": rho_ref,
            "rho_error_percent": rho_error,
            "cp_PR_J_kgK": cp_pr, "cp_HEOS_J_kgK": cp_ref,
            "cp_error_percent": cp_error,
            "delta_h_PR_from_anchor_J_kg": delta_h_pr,
            "delta_h_HEOS_from_anchor_J_kg": delta_h_ref,
            "delta_h_error_J_kg": dh_error,
            "delta_s_PR_from_anchor_J_kgK": delta_s_pr,
            "delta_s_HEOS_from_anchor_J_kgK": delta_s_ref,
            "delta_s_error_J_kgK": ds_error,
        }
        output_rows.append(comparison)
        by_group[group].append(comparison)

    if len(output_rows) != len(pr_rows):
        raise RuntimeError(
            f"PR property export contains {len(pr_rows) - len(output_rows)} failed states; "
            "inspect the C++ export before interpreting this comparison"
        )

    with args.output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    def mean_abs(values: list[float]) -> float:
        return sum(abs(value) for value in values) / len(values)

    print(f"CoolProp version: {CoolProp.__version__}; backend: HEOS::CO2")
    print(f"States compared: {len(output_rows)}; accuracy gate: none")
    for group, rows in by_group.items():
        if not rows:
            continue
        rho_errors = [float(row["rho_error_percent"]) for row in rows]
        cp_errors = [float(row["cp_error_percent"]) for row in rows]
        delta_h_errors = [float(row["delta_h_error_J_kg"]) for row in rows]
        delta_s_errors = [float(row["delta_s_error_J_kgK"]) for row in rows]
        print(f"{group}: n={len(rows)}")
        print(f"  density: MAPE={mean_abs(rho_errors):.3f}%, max_abs={max(map(abs, rho_errors)):.3f}%")
        print(f"  Cp: MAPE={mean_abs(cp_errors):.3f}%, max_abs={max(map(abs, cp_errors)):.3f}%")
        print(f"  datum-aligned enthalpy MAE={mean_abs(delta_h_errors):.2f} J/kg")
        print(f"  datum-aligned entropy MAE={mean_abs(delta_s_errors):.5f} J/(kg K)")
    print(f"Pointwise comparison: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
