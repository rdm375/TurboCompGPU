#pragma once

#include <turbocompgpu/vaneless_diffuser.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <utility>
#include <stdexcept>

namespace turbocompgpu {

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

namespace frictional_diffuser_integrator_detail {

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

} // namespace frictional_diffuser_integrator_detail

template<class GeometryProvider, class DerivativeProvider>
inline FrictionalDiffuserResult
integrate_frictional_vaneless_diffuser_profile(
    const thermomodel::Model& model,
    double inlet_radius,
    double outlet_radius,
    GeometryProvider&& geometry_provider,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options,
    DerivativeProvider&& derivative_provider)
{

    if (!(std::isfinite(inlet_radius) &&
          std::isfinite(outlet_radius) &&
          inlet_radius > 0.0 &&
          outlet_radius > inlet_radius &&
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

    const auto inlet_geometry = geometry_provider(inlet_radius);
    const double inlet_width = inlet_geometry.first;
    if (!(std::isfinite(inlet_width) && inlet_width > 0.0 &&
          std::isfinite(inlet_geometry.second))) {
        throw std::invalid_argument(
            "Invalid diffuser inlet geometry");
    }
    const double area2 =
        2.0 * std::numbers::pi *
        inlet_radius * inlet_width;

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

    frictional_diffuser_integrator_detail::Vector y{
        state2.p,
        state2.T,
        inlet.tangential_velocity
    };

    double entropy_integral = 0.0;
    double angular_momentum_integral = 0.0;

    const double dr =
        (outlet_radius - inlet_radius) /
        options.radial_steps;

    for (int i = 0; i < options.radial_steps; ++i) {
        const double r =
            inlet_radius + i * dr;

        auto rhs = [&](double radius, const frictional_diffuser_integrator_detail::Vector& state) {
            const auto [width, width_derivative] =
                geometry_provider(radius);
            if (!(std::isfinite(width) && width > 0.0 &&
                  std::isfinite(width_derivative))) {
                throw std::invalid_argument(
                    "Invalid diffuser width profile");
            }
            return derivative_provider(
                radius,
                width,
                width_derivative,
                inlet.mass_flow,
                options.darcy_friction_factor,
                options.derivative_relative_step,
                options.singularity_tolerance,
                state);
        };

        const auto k1 = rhs(r, y);
        const auto k2 = rhs(
            r + 0.5 * dr,
            frictional_diffuser_integrator_detail::add(y, k1.state, 0.5 * dr));
        const auto k3 = rhs(
            r + 0.5 * dr,
            frictional_diffuser_integrator_detail::add(y, k2.state, 0.5 * dr));
        const auto k4 = rhs(
            r + dr,
            frictional_diffuser_integrator_detail::add(y, k3.state, dr));

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

    const auto outlet_geometry = geometry_provider(outlet_radius);
    const double outlet_width = outlet_geometry.first;
    if (!(std::isfinite(outlet_width) && outlet_width > 0.0 &&
          std::isfinite(outlet_geometry.second))) {
        throw std::invalid_argument(
            "Invalid diffuser outlet geometry");
    }
    const double area3 =
        2.0 * std::numbers::pi *
        outlet_radius * outlet_width;

    const double u3 =
        inlet.mass_flow / (state3.rho * area3);

    const double h03 =
        state3.h + 0.5 * (u3 * u3 + y[2] * y[2]);

    const double L2 =
        inlet_radius * inlet.tangential_velocity;

    const double L3 =
        outlet_radius * y[2];

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


template<class DerivativeProvider>
inline FrictionalDiffuserResult
integrate_frictional_vaneless_diffuser(
    const thermomodel::Model& model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options,
    DerivativeProvider&& derivative_provider)
{
    if (!(std::isfinite(geometry.outlet_width) &&
          geometry.outlet_width > 0.0 &&
          geometry.outlet_blockage == 0.0)) {
        throw std::invalid_argument(
            "Invalid constant-width diffuser geometry");
    }

    auto profile = [&](double) {
        return std::pair{
            geometry.outlet_width, 0.0
        };
    };

    auto adapted_derivatives =
        [&](double radius,
            double width,
            double /*width_derivative*/,
            double mass_flow,
            double friction_factor,
            double relative_step,
            double singularity_tolerance,
            const frictional_diffuser_integrator_detail::Vector& y)
    {
        return derivative_provider(
            radius, width, mass_flow,
            friction_factor, relative_step,
            singularity_tolerance, y);
    };

    return integrate_frictional_vaneless_diffuser_profile(
        model,
        geometry.inlet_radius,
        geometry.outlet_radius,
        profile,
        inlet,
        options,
        adapted_derivatives);
}

} // namespace turbocompgpu
