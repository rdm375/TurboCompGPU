#pragma once

#include <turbocompgpu/vaneless_diffuser_friction.hpp>
#include <turbocompgpu/vaneless_diffuser_friction_analytical.hpp>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

namespace turbocompgpu {

struct VariableWidthDiffuserGeometry {
    double inlet_radius;
    double outlet_radius;

    std::function<double(double)> width;
    std::function<double(double)> width_derivative;
};

namespace variable_width_diffuser_detail {

template<class DerivativeProvider>
inline FrictionalDiffuserResult integrate(
    const thermomodel::Model& model,
    const VariableWidthDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options,
    DerivativeProvider&& derivatives)
{
    if (!(std::isfinite(geometry.inlet_radius) &&
          std::isfinite(geometry.outlet_radius) &&
          geometry.inlet_radius > 0.0 &&
          geometry.outlet_radius > geometry.inlet_radius &&
          geometry.width &&
          geometry.width_derivative)) {
        throw std::invalid_argument(
            "Invalid variable-width diffuser geometry");
    }

    auto profile = [&](double radius) {
        return std::pair{
            geometry.width(radius),
            geometry.width_derivative(radius)
        };
    };

    return integrate_frictional_vaneless_diffuser_profile(
        model,
        geometry.inlet_radius,
        geometry.outlet_radius,
        profile,
        inlet,
        options,
        derivatives);
}

} // namespace variable_width_diffuser_detail

inline FrictionalDiffuserResult
evaluate_variable_width_vaneless_diffuser(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives_model,
    const VariableWidthDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options = {})
{
    auto derivatives =
        [&](double radius,
            double width,
            double width_derivative,
            double mass_flow,
            double friction_factor,
            double relative_step,
            double singularity_tolerance,
            const frictional_diffuser_integrator_detail::Vector& y)
    {
        return frictional_diffuser_analytical_detail::evaluate_rhs(
            model,
            derivatives_model,
            radius,
            width,
            mass_flow,
            friction_factor,
            relative_step,
            singularity_tolerance,
            y,
            width_derivative);
    };

    return variable_width_diffuser_detail::integrate(
        model, geometry, inlet, options, derivatives);
}

inline FrictionalDiffuserResult
evaluate_variable_width_vaneless_diffuser_finite_difference(
    const thermomodel::Model& model,
    const VariableWidthDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options = {})
{
    auto derivatives =
        [&](double radius,
            double width,
            double width_derivative,
            double mass_flow,
            double friction_factor,
            double relative_step,
            double singularity_tolerance,
            const frictional_diffuser_integrator_detail::Vector& y)
    {
        return frictional_diffuser_detail::evaluate_rhs(
            model,
            radius,
            width,
            mass_flow,
            friction_factor,
            relative_step,
            singularity_tolerance,
            y,
            width_derivative);
    };

    return variable_width_diffuser_detail::integrate(
        model, geometry, inlet, options, derivatives);
}

} // namespace turbocompgpu
