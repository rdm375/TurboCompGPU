#!/usr/bin/env python3
"""Independent arithmetic and convention screen for Eckardt A loss CSV."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def fanning_factor(reynolds: float) -> float:
    low = 2.67 * reynolds ** -0.5
    high = 0.0622 * reynolds ** -0.2
    if reynolds <= 2.0e5:
        return low
    if reynolds >= 3.0e5:
        return high
    w = (reynolds - 2.0e5) / 1.0e5
    return low * (1.0 - w) + high * w


def viscosity_sutherland(temperature_k: float) -> float:
    t_ref, mu_ref, s = 273.15, 1.716e-5, 110.4
    return mu_ref * (temperature_k / t_ref) ** 1.5 * (t_ref + s) / (temperature_k + s)


def expected_geometry(angle_deg: float) -> tuple[float, float]:
    pi = math.pi
    z, d2, b2 = 20.0, 0.400, 0.026
    d1s, d1h, lz = 0.280, 0.120, 0.130
    b1s, b1h = math.radians(63.0), math.radians(40.0)
    beta2 = math.radians(angle_deg)
    dh_exit = d2 * math.cos(beta2) / (z / pi + d2 * math.cos(beta2) / b2)
    beta1_mean_cos = 0.5 * (math.cos(b1s) + math.cos(b1h))
    dh_inlet = 0.5 * (d1s + d1h) * beta1_mean_cos / (
        z / pi + (d1s + d1h) / (d1s - d1h) * beta1_mean_cos
    )
    dh = dh_exit + dh_inlet
    path = pi / 8.0 * (d2 - ((d1s + d1h) / 2.0) - b2 + 2.0 * lz) * (
        2.0 / (0.5 * (math.cos(b1s) + math.cos(b1h)) + math.cos(beta2))
    )
    return dh, path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="?", default="benchmarks/results-eckardt-impeller-a/loss_term_ablation.csv")
    args = parser.parse_args()
    rows = list(csv.DictReader(Path(args.csv).open(newline="")))
    if not rows:
        raise SystemExit("No rows found")

    max_re_rel = max_skin_rel = max_dh_rel = max_path_rel = 0.0
    max_throat_mass_rel = max_throat_rothalpy_rel = max_wbar_rel = 0.0
    max_throat_entropy_abs = 0.0
    max_throat_mach = 0.0
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["scenario"]].append(row)
        angle = float(row["impeller_correlation_exit_angle_deg"])
        dh_expected, path_expected = expected_geometry(angle)
        if row["scenario"] == "rectangular_endpoint_hydraulic_diameter":
            pi = math.pi
            r1s, r1h, r1m, r2, b2, z = 0.140, 0.060, math.sqrt(0.5 * (0.140**2 + 0.060**2)), 0.200, 0.026, 20.0
            inlet_pitch = 2.0 * pi * r1m / z - 0.00211
            exit_blockage = z * 0.00108 / (2.0 * pi * r2)
            exit_pitch = 2.0 * pi * r2 / z * (1.0 - exit_blockage)
            dh1 = 2.0 * (r1s-r1h) * inlet_pitch / ((r1s-r1h)+inlet_pitch)
            dh2 = 2.0 * b2 * exit_pitch / (b2+exit_pitch)
            dh_expected = 0.5 * (dh1 + dh2)
        dh = float(row["impeller_hydraulic_diameter_m"])
        path = float(row["impeller_blade_path_length_m"])
        reynolds = float(row["impeller_reynolds"])
        w = float(row["impeller_mean_relative_velocity_m_s"])
        rho = float(row["impeller_mean_density_kg_m3"])
        temp = float(row["impeller_mean_temperature_k"])
        cf = float(row["impeller_skin_fanning_factor"])
        mu = viscosity_sutherland(temp)
        re_expected = rho * w * dh / mu
        cf_expected = fanning_factor(reynolds)
        skin_expected = 4.0 * cf_expected * (path / dh) * w * w
        skin_actual = float(row["impeller_skin_friction_loss_j_kg"])
        max_re_rel = max(max_re_rel, abs(re_expected - reynolds) / reynolds)
        max_skin_rel = max(max_skin_rel, abs(skin_expected - skin_actual) / max(1.0, abs(skin_actual)))
        max_dh_rel = max(max_dh_rel, abs(dh_expected - dh) / dh)
        max_path_rel = max(max_path_rel, abs(path_expected - path) / path)

        flow = float(row["mass_flow_kg_s"])
        wth = float(row["impeller_throat_relative_velocity_m_s"])
        rho_th = float(row["impeller_throat_density_kg_m3"])
        area_th = float(row["impeller_throat_area_m2"])
        continuity = rho_th * area_th * wth - flow
        h_th = float(row["impeller_throat_enthalpy_j_kg"])
        ut = float(row["impeller_throat_blade_speed_m_s"])
        rothalpy = float(row["impeller_throat_rothalpy_j_kg"])
        energy = h_th + 0.5 * wth**2 - 0.5 * ut**2 - rothalpy
        max_throat_mass_rel = max(max_throat_mass_rel, abs(continuity) / max(1.0, flow))
        max_throat_rothalpy_rel = max(max_throat_rothalpy_rel, abs(energy) / max(1.0, abs(rothalpy)))
        max_throat_entropy_abs = max(max_throat_entropy_abs,
            abs(float(row["impeller_throat_entropy_residual_j_kg_k"])))
        mach = float(row["impeller_throat_mach"])
        max_throat_mach = max(max_throat_mach, mach)
        if mach >= 1.0:
            raise SystemExit(f"Throat solution is not subsonic in {row['scenario']}")

        w1 = float(row["impeller_inlet_mean_relative_velocity_m_s"])
        w1s = float(row["impeller_inlet_shroud_relative_velocity_m_s"])
        w2 = float(row["impeller_exit_relative_velocity_m_s"])
        second_speed = w1s if row["scenario"] == "shroud_inlet_speed_proxy" else wth
        wbar_expected = max(math.sqrt(0.5 * (w1*w1 + w2*w2)),
                            math.sqrt(0.5 * (second_speed*second_speed + w2*w2)))
        wbar = float(row["impeller_mean_relative_velocity_m_s"])
        max_wbar_rel = max(max_wbar_rel, abs(wbar_expected-wbar)/max(1.0, abs(wbar)))

        if abs(cf_expected - cf) > 1.0e-12 * max(1.0, abs(cf)):
            raise SystemExit(f"Fanning correlation mismatch in {row['scenario']}")

    print(f"Rows checked: {len(rows)}")
    print(f"Max relative residuals: Re={max_re_rel:.3e}, Galvas Dh={max_dh_rel:.3e}, Jansen path={max_path_rel:.3e}, skin-loss arithmetic={max_skin_rel:.3e}")
    print(f"Operating-throat checks: continuity={max_throat_mass_rel:.3e}, rothalpy={max_throat_rothalpy_rel:.3e}, entropy_abs={max_throat_entropy_abs:.3e} J/(kg K), Wbar={max_wbar_rel:.3e}, max_Mach={max_throat_mach:.6f}")
    if max(max_re_rel, max_dh_rel, max_path_rel, max_skin_rel,
           max_throat_mass_rel, max_throat_rothalpy_rel, max_wbar_rel) > 2.0e-10 or max_throat_entropy_abs > 2.0e-8:
        raise SystemExit("Arithmetic or operating-throat audit failed")

    print("\nScenario screening metrics (approximate efficiency digitization):")
    if any(float(r["impeller_mixing_bstar"]) < 1.0 for r in rows):
        print("WARNING: bstar < 1 rows are outside Johnston-Dean's sudden-expansion model domain.")
    else:
        print("Mixing bstar=1 is the model boundary; bstar>1 cases are inside the stated domain.")
    print("scenario,angle_deg,mean_Dh_mm,mean_Lb_mm,mean_skin_kJ_kg,mean_mix_kJ_kg,eta_bias_pp,eta_MAE_pp")
    wanted = {
        "all_losses_enthalpy_equivalent",
        "shroud_inlet_speed_proxy",
        "correlation_exit_angle_60deg_screen",
        "mixing_chi_0.05_bstar_1.00",
        "mixing_chi_0.10_bstar_1.00",
        "mixing_chi_0.20_bstar_1.00",
        "mixing_chi_0.25_bstar_1.00",
        "mixing_chi_0.30_bstar_1.00",
        "mixing_chi_0.15_bstar_1.02",
        "mixing_chi_0.15_bstar_1.10",
    }
    for name in sorted(wanted):
        group = groups.get(name)
        if not group:
            continue
        avg = lambda key: sum(float(r[key]) for r in group) / len(group)
        errors = [float(r["eta_error_percentage_points"]) for r in group]
        print(f"{name},{avg('impeller_correlation_exit_angle_deg'):.1f},"
              f"{1000*avg('impeller_hydraulic_diameter_m'):.3f},"
              f"{1000*avg('impeller_blade_path_length_m'):.3f},"
              f"{avg('impeller_skin_friction_loss_j_kg')/1000:.3f},"
              f"{avg('impeller_mixing_loss_j_kg')/1000:.3f},"
              f"{sum(errors)/len(errors):.3f},"
              f"{sum(abs(e) for e in errors)/len(errors):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
