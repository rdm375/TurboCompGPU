#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct ImpellerGeometry {
    double exit_diameter;          // m
    double exit_blade_angle_deg;   // degrees

    // Representative inlet meanline diameter.
    // Zero retains the legacy zero-prewhirl configuration.
    double inlet_diameter = 0.0;   // m
};

struct ImpellerOperatingPoint {
    double rotational_speed_rpm;
    double exit_meridional_velocity; // m/s
    double slip_factor;

    // Optional inlet conditions.
    double inlet_meridional_velocity = 0.0; // m/s
    double inlet_tangential_velocity = 0.0; // m/s
};

struct ImpellerResult {
    // Existing exit quantities.
    double tip_speed;                 // m/s
    double ideal_exit_whirl;          // m/s
    double actual_exit_whirl;         // m/s
    double absolute_exit_speed;       // m/s
    double relative_exit_speed;       // m/s
    double relative_exit_angle_deg;   // degrees
    double euler_specific_work;       // J/kg

    // New inlet quantities.
    double inlet_blade_speed;         // m/s
    double inlet_absolute_speed;      // m/s
    double inlet_relative_speed;      // m/s
    double inlet_relative_angle_deg;  // degrees
    double inlet_euler_work;          // J/kg
    double exit_euler_work;           // J/kg
};

inline ImpellerResult evaluate_impeller(
    const ImpellerGeometry& geometry,
    const ImpellerOperatingPoint& operating_point)
{
    const double D2 = geometry.exit_diameter;
    const double D1 = geometry.inlet_diameter;
    const double beta_deg = geometry.exit_blade_angle_deg;

    const double N = operating_point.rotational_speed_rpm;

    const double Vm2 =
        operating_point.exit_meridional_velocity;

    const double sigma = operating_point.slip_factor;

    const double Vm1 =
        operating_point.inlet_meridional_velocity;

    const double Vtheta1 =
        operating_point.inlet_tangential_velocity;

    if (!std::isfinite(D2) || D2 <= 0.0) {
        throw std::invalid_argument(
            "Impeller exit diameter must be positive and finite");
    }

    if (!std::isfinite(D1) || D1 < 0.0 || D1 >= D2) {
        throw std::invalid_argument(
            "Impeller inlet diameter must be nonnegative and less than exit diameter");
    }

    if (!std::isfinite(beta_deg) ||
        beta_deg <= 0.0 ||
        beta_deg > 90.0) {
        throw std::invalid_argument(
            "Impeller exit blade angle must be in (0,90] degrees");
    }

    if (!std::isfinite(N) || N <= 0.0) {
        throw std::invalid_argument(
            "Rotational speed must be positive and finite");
    }

    if (!std::isfinite(Vm2) || Vm2 < 0.0) {
        throw std::invalid_argument(
            "Exit meridional velocity must be nonnegative and finite");
    }

    if (!std::isfinite(Vm1) || Vm1 < 0.0) {
        throw std::invalid_argument(
            "Inlet meridional velocity must be nonnegative and finite");
    }

    if (!std::isfinite(Vtheta1)) {
        throw std::invalid_argument(
            "Inlet tangential velocity must be finite");
    }

    if (D1 == 0.0 && Vtheta1 != 0.0) {
        throw std::invalid_argument(
            "Nonzero inlet prewhirl requires a positive inlet diameter");
    }

    if (!std::isfinite(sigma) ||
        sigma <= 0.0 ||
        sigma > 1.0) {
        throw std::invalid_argument(
            "Slip factor must be in (0,1]");
    }

    const double beta =
        beta_deg * std::numbers::pi / 180.0;

    const double U1 =
        std::numbers::pi * D1 * N / 60.0;

    const double U2 =
        std::numbers::pi * D2 * N / 60.0;

    const double cot_beta =
        std::cos(beta) / std::sin(beta);

    const double Vtheta2_ideal =
        U2 - Vm2 * cot_beta;

    const double Vtheta2 =
        sigma * Vtheta2_ideal;

    // Exit velocity triangle.
    const double V2 =
        std::hypot(Vm2, Vtheta2);

    const double Wtheta2 =
        Vtheta2 - U2;

    const double W2 =
        std::hypot(Vm2, Wtheta2);

    const double beta_flow2 =
        std::atan2(Vm2, -Wtheta2) *
        180.0 / std::numbers::pi;

    // Inlet velocity triangle.
    const double V1 =
        std::hypot(Vm1, Vtheta1);

    const double Wtheta1 =
        Vtheta1 - U1;

    const double W1 =
        std::hypot(Vm1, Wtheta1);

    const double beta_flow1 =
        std::atan2(Vm1, -Wtheta1) *
        180.0 / std::numbers::pi;

    // Complete Euler turbomachinery equation.
    const double inlet_work =
        U1 * Vtheta1;

    const double exit_work =
        U2 * Vtheta2;

    const double work =
        exit_work - inlet_work;

    if (!std::isfinite(U1) ||
        !std::isfinite(U2) ||
        !std::isfinite(Vtheta2_ideal) ||
        !std::isfinite(Vtheta2) ||
        !std::isfinite(V1) ||
        !std::isfinite(W1) ||
        !std::isfinite(V2) ||
        !std::isfinite(W2) ||
        !std::isfinite(beta_flow1) ||
        !std::isfinite(beta_flow2) ||
        !std::isfinite(inlet_work) ||
        !std::isfinite(exit_work) ||
        !std::isfinite(work)) {
        throw std::runtime_error(
            "Non-finite impeller aerodynamic result");
    }

    return {
        U2,
        Vtheta2_ideal,
        Vtheta2,
        V2,
        W2,
        beta_flow2,
        work,
        U1,
        V1,
        W1,
        beta_flow1,
        inlet_work,
        exit_work
    };
}

} // namespace turbocompgpu
