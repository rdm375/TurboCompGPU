#pragma once

#include <turbocompgpu/operating_incidence.hpp>
#include <turbocompgpu/operating_point.hpp>
#include <turbocompgpu/station.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

struct IncidenceOperatingPointResult {
    OperatingPointResult solution;

    OperatingIncidenceResult incidence;

    double inlet_stagnation_temperature;
    double baseline_entropy_generation;
    double incidence_entropy_generation;
    double total_entropy_generation;
};

inline IncidenceOperatingPointResult
evaluate_incidence_operating_point(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const IncidenceLossParameters& incidence_parameters,
    const NonlinearSolverOptions& options = {},
    double validation_tolerance = 1.0e-8)
{
    const auto& operating = request.operating;

    if (!std::isfinite(operating.entropy_generation) ||
        operating.entropy_generation < 0.0) {
        throw std::invalid_argument(
            "Invalid baseline entropy generation");
    }

    // Evaluate the inlet state once.
    const auto inlet_evaluation = model.evaluate_Tp(
        operating.inlet_temperature,
        operating.inlet_pressure
    );

    if (!inlet_evaluation) {
        throw std::runtime_error(
            "Inlet thermodynamic evaluation failed: " +
            inlet_evaluation.message
        );
    }

    const auto& inlet = inlet_evaluation.value;

    // The incidence loss depends only on inlet
    // conditions, not the nonlinear outlet unknowns.
    const OperatingIncidenceConditions conditions {
        operating.mass_flow,
        operating.rotational_speed_rpm,
        inlet.rho,
        operating.inlet_prewhirl
    };

    const auto incidence =
        evaluate_operating_incidence_loss(
            geometry,
            conditions,
            incidence_parameters
        );

    // Use the existing thermodynamic stagnation-state
    // evaluator, preserving EOS consistency.
    const FlowStation inlet_station {
        inlet,
        {
            incidence.inlet_meridional_velocity,
            operating.inlet_prewhirl
        }
    };

    const auto stagnation =
        evaluate_stagnation_state(
            model,
            inlet_station
        );

    const double Tt1 =
        stagnation.stagnation_state.T;

    if (!std::isfinite(Tt1) || Tt1 <= 0.0) {
        throw std::runtime_error(
            "Invalid inlet stagnation temperature"
        );
    }

    // First-order engineering approximation.
    const double incidence_entropy =
        incidence.loss.specific_loss / Tt1;

    const double total_entropy =
        operating.entropy_generation +
        incidence_entropy;

    if (!std::isfinite(incidence_entropy) ||
        incidence_entropy < 0.0 ||
        !std::isfinite(total_entropy)) {
        throw std::runtime_error(
            "Invalid incidence entropy contribution"
        );
    }

    // Preserve the existing nonlinear solver API.
    OperatingPointRequest adjusted = request;

    adjusted.operating.entropy_generation =
        total_entropy;

    const auto solution = evaluate_operating_point(
        model,
        geometry,
        adjusted,
        scales,
        options,
        validation_tolerance
    );

    return {
        solution,
        incidence,
        Tt1,
        operating.entropy_generation,
        incidence_entropy,
        total_entropy
    };
}

} // namespace turbocompgpu
