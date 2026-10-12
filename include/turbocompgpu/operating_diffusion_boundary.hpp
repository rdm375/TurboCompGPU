#pragma once

#include <turbocompgpu/diffusion_boundary.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>
#include <turbocompgpu/fast_nonlinear_residual.hpp>
#include <turbocompgpu/incidence_operating_point.hpp>
#include <turbocompgpu/operating_diffusion.hpp>

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct OperatingDiffusionBoundaryEvaluation {
    double diffusion_ratio;
    double specific_loss;
    double inlet_meridional_velocity;
    double outlet_meridional_velocity;
    double outlet_whirl_velocity;

    IncidenceOperatingPointResult incidence_point;
};

inline OperatingDiffusionBoundaryEvaluation
evaluate_operating_diffusion_boundary(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const IncidenceLossParameters& incidence_parameters,
    const DiffusionLossParameters& diffusion_parameters,
    const NonlinearSolverOptions& options = {})
{
    const auto incidence = evaluate_incidence_operating_point(
        model, geometry, request, scales,
        incidence_parameters, options
    );

    const auto& solution = incidence.solution;

    if (!solution.converged ||
        !solution.independently_validated ||
        !solution.basic_physics_valid) {
        throw std::runtime_error(
            "Invalid incidence operating point in boundary evaluation"
        );
    }

    // The fast residual must use the entropy generation
    // of the converged incidence-only solution.
    const auto context = make_fast_residual_context(
        model, geometry, solution.operating, scales
    );

    const NonlinearUnknowns unknowns {
        solution.outlet_temperature,
        std::log(solution.outlet_pressure / scales.pressure)
    };

    const auto flow = evaluate_fast_nonlinear_residual(
        context, unknowns
    );

    if (!std::isfinite(flow.residual_norm) ||
        flow.residual_norm > 1.0e-8) {
        throw std::runtime_error(
            "Fast residual reconstruction failed"
        );
    }

    const OperatingDiffusionConditions conditions {
        incidence.incidence.inlet_meridional_velocity,
        flow.outlet_meridional_velocity,
        request.operating.inlet_prewhirl,
        flow.outlet_whirl,
        request.operating.rotational_speed_rpm
    };

    const auto diffusion = evaluate_operating_diffusion_loss(
        geometry, conditions, diffusion_parameters
    );

    return {
        diffusion.loss.diffusion_ratio,
        diffusion.loss.specific_loss,
        incidence.incidence.inlet_meridional_velocity,
        flow.outlet_meridional_velocity,
        flow.outlet_whirl,
        incidence
    };
}

} // namespace turbocompgpu
