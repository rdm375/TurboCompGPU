#pragma once

#include <turbocompgpu/vaneless_diffuser.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

// Stage 8C.1b: steady, adiabatic, frictional vaneless diffuser.
//
// Initial restrictions:
//   * constant width
//   * zero blockage
//   * constant Darcy friction factor
//   * outward radial flow
//   * smooth, nonsingular solution
//
// Unknowns: static pressure, temperature, tangential velocity.
// Radial velocity follows from continuity.
//
// This is a CPU reference integrator, not a choking model.

struct FrictionalDiffuserOptions {
    double darcy_friction_factor = 0.0;
    int radial_steps = 512;
    double derivative_relative_step = 1.0e-5;
    double singularity_tolerance = 1.0e-10;
};

struct FrictionalDiffuserResult {
    thermomodel::State inlet_static;
    thermomodel::State outlet_static;

    double inlet_radial_velocity = 0.0;
    double inlet_tangential_velocity = 0.0;

    double outlet_radial_velocity = 0.0;
    double outlet_tangential_velocity = 0.0;

    double mass_flow = 0.0;

    double mass_residual = 0.0;
    double stagnation_enthalpy_residual = 0.0;

    double entropy_change = 0.0;
    double integrated_entropy_production = 0.0;
    double entropy_balance_residual = 0.0;

    double angular_momentum_change = 0.0;
    double integrated_wall_torque_per_mass = 0.0;
    double angular_momentum_balance_residual = 0.0;

    int radial_steps = 0;
};

namespace frictional_diffuser_detail {

using Vector = std::array<double, 3>;

inline Vector add(
    const Vector& a,
    const Vector& b,
    double factor)
{
    return {
        a[0] + factor * b[0],
        a[1] + factor * b[1],
        a[2] + factor * b[2]
    };
}

struct Derivatives {
    Vector state;
    double entropy;
    double angular_momentum;
};

inline Derivatives evaluate_rhs(
    const thermomodel::Model& model,
    double radius,
    double width,
    double mass_flow,
    double friction_factor,
    double relative_step,
    double singularity_tolerance,
    const Vector& y)
{
    const double p = y[0];
    const double T = y[1];
    const double v = y[2];

    if (!(std::isfinite(p) && std::isfinite(T) &&
          std::isfinite(v) && p > 0.0 && T > 0.0)) {
        throw std::runtime_error(
            "Frictional diffuser left thermodynamic domain");
    }

    const auto state = evaluate(model, T, p);

    if (!(std::isfinite(state.rho) && state.rho > 0.0 &&
          std::isfinite(state.h))) {
        throw std::runtime_error(
            "Invalid frictional diffuser thermodynamic state");
    }

    const double area =
        2.0 * std::numbers::pi * radius * width;

    const double u = mass_flow / (state.rho * area);
    const double V = std::hypot(u, v);

    // Central differences in the primitive thermodynamic variables.
    const double dp = relative_step * p;
    const double dT = relative_step * T;

    if (!(dp > 0.0 && dT > 0.0 && p > dp && T > dT)) {
        throw std::runtime_error(
            "Invalid thermodynamic differentiation interval");
    }

    const auto pp = evaluate(model, T, p + dp);
    const auto pm = evaluate(model, T, p - dp);
    const auto Tp = evaluate(model, T + dT, p);
    const auto Tm = evaluate(model, T - dT, p);

    const double rho_p =
        (pp.rho - pm.rho) / (2.0 * dp);

    const double rho_T =
        (Tp.rho - Tm.rho) / (2.0 * dT);

    const double h_p =
        (pp.h - pm.h) / (2.0 * dp);

    const double h_T =
        (Tp.h - Tm.h) / (2.0 * dT);

    // Wall shear: tau = f_D rho V^2 / 8.
    // Two annular walls give a friction acceleration
    // coefficient f_D V / (4 b).
    const double friction =
        friction_factor * V / (4.0 * width);

    const double dv =
        -v / radius - friction * v / u;

    // Radial momentum and stagnation-enthalpy conservation
    // form a 2x2 system for dp/dr and dT/dr.
    const double A11 =
        (1.0 - u * u * rho_p) / state.rho;

    const double A12 =
        -u * u * rho_T / state.rho;

    const double A21 =
        h_p - u * u * rho_p / state.rho;

    const double A22 =
        h_T - u * u * rho_T / state.rho;

    const double B1 =
        (u * u + v * v) / radius - friction * u;

    const double B2 =
        u * u / radius - v * dv;

    const double determinant =
        A11 * A22 - A12 * A21;

    const double determinant_scale =
        std::abs(A11 * A22) +
        std::abs(A12 * A21);

    if (!(std::isfinite(determinant) &&
          std::isfinite(determinant_scale) &&
          determinant_scale > 0.0 &&
          std::abs(determinant) >
              singularity_tolerance * determinant_scale)) {
        throw std::runtime_error(
            "Frictional diffuser derivative system is singular");
    }

    const double dpdr =
        (B1 * A22 - A12 * B2) / determinant;

    const double dTdr =
        (A11 * B2 - B1 * A21) / determinant;

    // Independent integral diagnostics.
    const double dsdr =
        friction_factor * V * V * V /
        (4.0 * width * u * T);

    const double dLdr =
        -friction * radius * v / u;

    if (!(std::isfinite(dpdr) &&
          std::isfinite(dTdr) &&
          std::isfinite(dv) &&
          std::isfinite(dsdr) &&
          std::isfinite(dLdr))) {
        throw std::runtime_error(
            "Non-finite frictional diffuser derivative");
    }

    return {{dpdr, dTdr, dv}, dsdr, dLdr};
}

} // namespace frictional_diffuser_detail

inline FrictionalDiffuserResult evaluate_frictional_vaneless_diffuser(
    const thermomodel::Model& model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options = {})
{
    using namespace frictional_diffuser_detail;

    if (!(std::isfinite(geometry.inlet_radius) &&
          std::isfinite(geometry.outlet_radius) &&
          std::isfinite(geometry.outlet_width) &&
          geometry.inlet_radius > 0.0 &&
          geometry.outlet_radius > geometry.inlet_radius &&
          geometry.outlet_width > 0.0 &&
          geometry.outlet_blockage == 0.0 &&
          std::isfinite(inlet.mass_flow) &&
          inlet.mass_flow > 0.0 &&
          std::isfinite(inlet.meridional_velocity) &&
          inlet.meridional_velocity > 0.0 &&
          std::isfinite(inlet.tangential_velocity) &&
          std::isfinite(options.darcy_friction_factor) &&
          options.darcy_friction_factor >= 0.0 &&
          options.radial_steps > 0 &&
          options.radial_steps <= 1000000 &&
          std::isfinite(options.derivative_relative_step) &&
          options.derivative_relative_step > 0.0 &&
          options.derivative_relative_step < 0.1 &&
          std::isfinite(options.singularity_tolerance) &&
          options.singularity_tolerance > 0.0)) {
        throw std::invalid_argument(
            "Invalid frictional diffuser input or options");
    }

    const auto state2 = evaluate(
        model, inlet.static_state.T, inlet.static_state.p);

    const double area2 =
        2.0 * std::numbers::pi *
        geometry.inlet_radius * geometry.outlet_width;

    const double expected_mass =
        state2.rho * area2 * inlet.meridional_velocity;

    if (std::abs(expected_mass - inlet.mass_flow) >
        1.0e-8 * std::max(1.0, inlet.mass_flow)) {
        throw std::invalid_argument(
            "Diffuser inlet violates continuity");
    }

    const double h02 =
        state2.h +
        0.5 * (
            inlet.meridional_velocity *
                inlet.meridional_velocity +
            inlet.tangential_velocity *
                inlet.tangential_velocity);

    Vector y{
        state2.p,
        state2.T,
        inlet.tangential_velocity
    };

    double entropy_integral = 0.0;
    double angular_momentum_integral = 0.0;

    const double dr =
        (geometry.outlet_radius - geometry.inlet_radius) /
        options.radial_steps;

    for (int i = 0; i < options.radial_steps; ++i) {
        const double r =
            geometry.inlet_radius + i * dr;

        auto rhs = [&](double radius, const Vector& state) {
            return evaluate_rhs(
                model,
                radius,
                geometry.outlet_width,
                inlet.mass_flow,
                options.darcy_friction_factor,
                options.derivative_relative_step,
                options.singularity_tolerance,
                state);
        };

        const auto k1 = rhs(r, y);
        const auto k2 = rhs(
            r + 0.5 * dr,
            add(y, k1.state, 0.5 * dr));
        const auto k3 = rhs(
            r + 0.5 * dr,
            add(y, k2.state, 0.5 * dr));
        const auto k4 = rhs(
            r + dr,
            add(y, k3.state, dr));

        for (int j = 0; j < 3; ++j) {
            y[j] += dr / 6.0 * (
                k1.state[j] +
                2.0 * k2.state[j] +
                2.0 * k3.state[j] +
                k4.state[j]);
        }

        entropy_integral += dr / 6.0 * (
            k1.entropy +
            2.0 * k2.entropy +
            2.0 * k3.entropy +
            k4.entropy);

        angular_momentum_integral += dr / 6.0 * (
            k1.angular_momentum +
            2.0 * k2.angular_momentum +
            2.0 * k3.angular_momentum +
            k4.angular_momentum);
    }

    const auto state3 = evaluate(model, y[1], y[0]);

    const double area3 =
        2.0 * std::numbers::pi *
        geometry.outlet_radius * geometry.outlet_width;

    const double u3 =
        inlet.mass_flow / (state3.rho * area3);

    const double h03 =
        state3.h + 0.5 * (u3 * u3 + y[2] * y[2]);

    const double L2 =
        geometry.inlet_radius * inlet.tangential_velocity;

    const double L3 =
        geometry.outlet_radius * y[2];

    FrictionalDiffuserResult result{};

    result.inlet_static = state2;
    result.outlet_static = state3;

    result.inlet_radial_velocity =
        inlet.meridional_velocity;
    result.inlet_tangential_velocity =
        inlet.tangential_velocity;

    result.outlet_radial_velocity = u3;
    result.outlet_tangential_velocity = y[2];

    result.mass_flow = inlet.mass_flow;

    result.mass_residual =
        state3.rho * area3 * u3 - inlet.mass_flow;

    result.stagnation_enthalpy_residual = h03 - h02;

    result.entropy_change = state3.s - state2.s;
    result.integrated_entropy_production = entropy_integral;

    result.entropy_balance_residual =
        result.entropy_change - entropy_integral;

    result.angular_momentum_change = L3 - L2;

    result.integrated_wall_torque_per_mass =
        angular_momentum_integral;

    result.angular_momentum_balance_residual =
        result.angular_momentum_change -
        angular_momentum_integral;

    result.radial_steps = options.radial_steps;

    return result;
}

} // namespace turbocompgpu
