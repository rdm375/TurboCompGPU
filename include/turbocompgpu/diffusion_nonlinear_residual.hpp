#pragma once

#include <turbocompgpu/fast_nonlinear_residual.hpp>
#include <turbocompgpu/operating_diffusion.hpp>

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct CoupledDiffusionResidual {
    FastNonlinearResidualResult flow;
    OperatingDiffusionResult diffusion;

    double fixed_entropy_generation;
    double diffusion_entropy_generation;
    double total_entropy_generation;

    // Required by the existing generic Newton solver.
    double normalized_energy;
    double normalized_entropy;
    double residual_norm;
};

inline CoupledDiffusionResidual
evaluate_coupled_diffusion_residual(
    const FastResidualContext& context,
    const NonlinearUnknowns& unknowns,
    double inlet_stagnation_temperature,
    const DiffusionLossParameters& parameters)
{
    if (!std::isfinite(inlet_stagnation_temperature) ||
        inlet_stagnation_temperature <= 0.0) {
        throw std::invalid_argument(
            "Invalid inlet stagnation temperature");
    }

    // The existing evaluator computes the outlet EOS state,
    // continuity velocity, slip-predicted whirl and Euler work.
    auto flow = evaluate_fast_nonlinear_residual(
        context, unknowns);

    const OperatingDiffusionConditions conditions {
        context.inlet_meridional_velocity,
        flow.outlet_meridional_velocity,
        context.operating.inlet_prewhirl,
        flow.outlet_whirl,
        context.operating.rotational_speed_rpm
    };

    const auto diffusion =
        evaluate_operating_diffusion_loss(
            context.geometry, conditions, parameters);

    const double diffusion_entropy =
        diffusion.loss.specific_loss /
        inlet_stagnation_temperature;

    const double fixed_entropy =
        context.operating.entropy_generation;

    const double total_entropy =
        fixed_entropy + diffusion_entropy;

    if (!std::isfinite(diffusion_entropy) ||
        !std::isfinite(total_entropy)) {
        throw std::runtime_error(
            "Non-finite coupled diffusion entropy");
    }

    // The original residual already subtracts fixed_entropy.
    // Only the state-dependent contribution is added here.
    flow.entropy_residual -= diffusion_entropy;

    flow.normalized_entropy =
        flow.entropy_residual /
        context.scales.specific_entropy;

    flow.residual_norm = std::hypot(
        flow.normalized_energy,
        flow.normalized_entropy);

    if (!std::isfinite(flow.residual_norm)) {
        throw std::runtime_error(
            "Non-finite coupled diffusion residual");
    }

    return {
        flow,
        diffusion,
        fixed_entropy,
        diffusion_entropy,
        total_entropy,
        flow.normalized_energy,
        flow.normalized_entropy,
        flow.residual_norm
    };
}

} // namespace turbocompgpu
