#pragma once

#include <turbocompgpu/vaneless_diffuser_friction_integrator.hpp>
#include <thermomodel/derivatives.hpp>
#include <friction/colebrook.hpp>
#include <turbocompgpu/iso20765_5_viscosity.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

namespace frictional_diffuser_analytical_detail {

inline frictional_diffuser_integrator_detail::Derivatives evaluate_rhs(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives_model,
    double radius,
    double width,
    double mass_flow,
    double friction_factor,
    double relative_step,
    double singularity_tolerance,
    const frictional_diffuser_integrator_detail::Vector& y,
    double width_derivative = 0.0)
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

    // Analytical EOS derivatives.
    const auto derivatives = derivatives_model.derivatives_Tp(T, p);

    if (!derivatives) {
        throw std::runtime_error(
            "Analytical diffuser EOS derivatives failed: " +
            derivatives.message);
    }

    const auto& d = derivatives.value;

    const double rho_p = d.drho_dp_T;
    const double rho_T = d.drho_dT_p;

    // Gibbs relation: dh = T ds + dp/rho.
    const double h_p =
        1.0 / state.rho + T * d.ds_dp_T;

    const double h_T =
        T * d.ds_dT_p;

    if (!(std::isfinite(rho_p) &&
          std::isfinite(rho_T) &&
          std::isfinite(h_p) &&
          std::isfinite(h_T) &&
          h_T > 0.0)) {
        throw std::runtime_error(
            "Invalid analytical diffuser derivatives");
    }

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

    // Logarithmic flow-area gradient:
    // A(r) = 2*pi*r*b(r).
    const double area_gradient =
        1.0 / radius + width_derivative / width;

    const double B1 =
        u * u * area_gradient +
        v * v / radius - friction * u;

    const double B2 =
        u * u * area_gradient - v * dv;

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

} // namespace frictional_diffuser_analytical_detail

inline FrictionalDiffuserResult
evaluate_frictional_vaneless_diffuser_analytical(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives_model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const FrictionalDiffuserOptions& options = {})
{
    auto provider = [&](double radius,
                        double width,
                        double mass_flow,
                        double friction_factor,
                        double relative_step,
                        double singularity_tolerance,
                        const frictional_diffuser_integrator_detail::Vector& y)
    {
        return frictional_diffuser_analytical_detail::evaluate_rhs(
            model, derivatives_model, radius, width,
            mass_flow, friction_factor, relative_step,
            singularity_tolerance, y);
    };

    return integrate_frictional_vaneless_diffuser(
        model, geometry, inlet, options, provider);
}

// Local Colebrook closure for the analytical diffuser integration. The local
// PR density and velocity are updated at each RK4 stage; viscosity is either
// evaluated from ISO 20765-5 at the local state or supplied as a fixed fallback.
// For a parallel-plate annulus, Dh = 2 b. This CPU reference mirrors the GPU
// stage-local closure.
inline FrictionalDiffuserResult
evaluate_frictional_vaneless_diffuser_local_colebrook_analytical(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives_model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    const transport::Iso20765Parameters& viscosity_model,
    bool use_iso20765_viscosity,
    double fixed_dynamic_viscosity_Pa_s,
    double relative_roughness,
    const FrictionalDiffuserOptions& options = {})
{
    if ((!use_iso20765_viscosity &&
          !(std::isfinite(fixed_dynamic_viscosity_Pa_s) &&
            fixed_dynamic_viscosity_Pa_s > 0.0)) ||
          !std::isfinite(relative_roughness) || relative_roughness < 0.0 ||
          relative_roughness >= 3.7) {
        throw std::invalid_argument("Invalid local Colebrook inputs");
    }

    auto provider = [&](double radius,
                        double width,
                        double mass_flow,
                        double,
                        double relative_step,
                        double singularity_tolerance,
                        const frictional_diffuser_integrator_detail::Vector& y)
    {
        const auto state = evaluate(model, y[1], y[0]);
        const double area = 2.0 * std::numbers::pi * radius * width;
        const double u = mass_flow / (state.rho * area);
        const double speed = std::hypot(u, y[2]);
        const double dynamic_viscosity = use_iso20765_viscosity
            ? transport::iso20765_5_dynamic_viscosity_Pa_s(
                  viscosity_model, y[1], y[0], state.rho)
            : fixed_dynamic_viscosity_Pa_s;
        if (!(std::isfinite(dynamic_viscosity) && dynamic_viscosity > 0.0))
            throw std::runtime_error("Local ISO 20765 viscosity evaluation failed");
        const double reynolds = state.rho * speed * (2.0 * width) /
            dynamic_viscosity;
        const auto friction = friction::colebrook_turbulent(
            reynolds, relative_roughness);
        if (!(std::isfinite(friction.darcy_friction_factor) &&
              friction.darcy_friction_factor > 0.0)) {
            throw std::runtime_error("Local CPU Colebrook evaluation failed");
        }
        return frictional_diffuser_analytical_detail::evaluate_rhs(
            model, derivatives_model, radius, width, mass_flow,
            friction.darcy_friction_factor, relative_step,
            singularity_tolerance, y, 0.0);
    };

    FrictionalDiffuserOptions integration_options = options;
    integration_options.darcy_friction_factor = 0.0;
    return integrate_frictional_vaneless_diffuser(
        model, geometry, inlet, integration_options, provider);
}

} // namespace turbocompgpu
