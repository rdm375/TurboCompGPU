#pragma once

#include <turbocompgpu/diffusion_nonlinear_residual.hpp>
#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/incidence_operating_point.hpp>
#include <turbocompgpu/nonlinear_solver.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

using CoupledDiffusionSolverResult =
    BasicNonlinearSolverResult<CoupledDiffusionResidual>;

struct DiffusionOperatingPointResult {
    CoupledDiffusionSolverResult solver;

    OperatingPointResult validation;

    OperatingIncidenceResult incidence;
    OperatingDiffusionResult diffusion;

    double inlet_stagnation_temperature;
    double baseline_entropy_generation;
    double incidence_entropy_generation;
    double diffusion_entropy_generation;
    double total_entropy_generation;

    double coupled_residual_norm;
    bool independently_validated;
};

inline DiffusionOperatingPointResult
evaluate_diffusion_operating_point(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const IncidenceLossParameters& incidence_parameters,
    const DiffusionLossParameters& diffusion_parameters,
    const NonlinearSolverOptions& options = {},
    double validation_tolerance = 1.0e-8)
{
    if (!std::isfinite(validation_tolerance) ||
        validation_tolerance <= 0.0) {
        throw std::invalid_argument(
            "Invalid diffusion validation tolerance");
    }

    const auto& operating = request.operating;

    if (!std::isfinite(operating.entropy_generation) ||
        operating.entropy_generation < 0.0) {
        throw std::invalid_argument(
            "Invalid baseline entropy generation");
    }

    // Inlet thermodynamics and incidence are independent
    // of the nonlinear outlet unknowns.
    const auto inlet_evaluation = model.evaluate_Tp(
        operating.inlet_temperature,
        operating.inlet_pressure);

    if (!inlet_evaluation) {
        throw std::runtime_error(
            "Diffusion inlet evaluation failed: " +
            inlet_evaluation.message);
    }

    const OperatingIncidenceConditions incidence_conditions {
        operating.mass_flow,
        operating.rotational_speed_rpm,
        inlet_evaluation.value.rho,
        operating.inlet_prewhirl
    };

    const auto incidence =
        evaluate_operating_incidence_loss(
            geometry,
            incidence_conditions,
            incidence_parameters);

    const FlowStation inlet_station {
        inlet_evaluation.value,
        {
            incidence.inlet_meridional_velocity,
            operating.inlet_prewhirl
        }
    };

    const auto stagnation =
        evaluate_stagnation_state(model, inlet_station);

    const double Tt1 =
        stagnation.stagnation_state.T;

    if (!std::isfinite(Tt1) || Tt1 <= 0.0) {
        throw std::runtime_error(
            "Invalid inlet stagnation temperature");
    }

    const double incidence_entropy =
        incidence.loss.specific_loss / Tt1;

    OperatingPointRequest adjusted = request;

    adjusted.operating.entropy_generation +=
        incidence_entropy;

    // The fast residual context retains the fixed entropy
    // contribution. Diffusion is added at every trial state.
    const auto context = make_fast_residual_context(
        model,
        geometry,
        adjusted.operating,
        scales);

    const auto evaluator =
        [&](const NonlinearUnknowns& unknowns) {
            return evaluate_coupled_diffusion_residual(
                context,
                unknowns,
                Tt1,
                diffusion_parameters);
        };

    const auto solved =
        solve_nonlinear_compressor_with_evaluator(
            evaluator,
            adjusted.initial_guess,
            options);

    if (!solved.converged() || !solved.solution) {
        throw std::runtime_error(
            "Coupled diffusion solver failed: " +
            solved.message);
    }

    const auto& coupled = *solved.solution;

    const double total_entropy =
        coupled.total_entropy_generation;

    // Independent reconstruction using the original
    // complete physical-system residual evaluator.
    OperatingPointRequest final_request = adjusted;

    final_request.operating.entropy_generation =
        total_entropy;

    final_request.initial_guess =
        coupled.flow.unknowns;

    // The final state is already converged. Re-solving here
    // is unnecessary; evaluate the full residual directly.
    const auto full = evaluate_nonlinear_residual(
        model,
        geometry,
        final_request.operating,
        coupled.flow.unknowns,
        scales);

    OperatingPointResult validation;

    validation.operating =
        final_request.operating;

    validation.status =
        NonlinearSolveStatus::Converged;

    validation.converged = true;

    validation.outlet_temperature =
        coupled.flow.unknowns.outlet_temperature;

    validation.outlet_pressure =
        coupled.flow.outlet_pressure;

    validation.pressure_ratio =
        validation.outlet_pressure /
        operating.inlet_pressure;

    validation.specific_work =
        coupled.flow.euler_specific_work;

    validation.shaft_power =
        operating.mass_flow *
        validation.specific_work;

    validation.torque =
        coupled.flow.shaft_torque;

    validation.entropy_generation =
        total_entropy;

    validation.inlet_density =
        context.inlet_state.rho;

    const auto outlet_evaluation = model.evaluate_Tp(
        validation.outlet_temperature,
        validation.outlet_pressure);

    if (!outlet_evaluation) {
        throw std::runtime_error(
            "Final diffusion outlet evaluation failed");
    }

    validation.outlet_density =
        outlet_evaluation.value.rho;

    validation.residual_norm =
        coupled.residual_norm;

    validation.validation_residual_norm =
        full.residual_norm;

    validation.iterations =
        solved.iterations;

    validation.evaluations =
        solved.residual_evaluations;

    validation.jacobians =
        solved.jacobian_evaluations;

    validation.rejected_steps =
        solved.rejected_steps;

    const double omega =
        context.angular_speed;

    const double torque_power =
        validation.torque * omega;

    const double power_scale =
        std::max(1.0, std::abs(validation.shaft_power));

    const bool power_consistent =
        std::abs(
            validation.shaft_power - torque_power
        ) <= validation_tolerance * power_scale;

    const bool independent =
        std::isfinite(full.residual_norm) &&
        full.residual_norm <= validation_tolerance &&
        std::abs(
            full.euler_specific_work -
            coupled.flow.euler_specific_work
        ) <= validation_tolerance *
            std::max(
                1.0,
                std::abs(full.euler_specific_work)
            ) &&
        std::abs(
            full.shaft_torque -
            coupled.flow.shaft_torque
        ) <= validation_tolerance *
            std::max(
                1.0,
                std::abs(full.shaft_torque)
            );

    validation.independently_validated =
        independent;

    validation.basic_physics_valid =
        independent &&
        power_consistent &&
        std::isfinite(validation.outlet_temperature) &&
        validation.outlet_temperature > 0.0 &&
        std::isfinite(validation.outlet_pressure) &&
        validation.outlet_pressure > 0.0 &&
        std::isfinite(validation.outlet_density) &&
        validation.outlet_density > 0.0 &&
        validation.pressure_ratio > 1.0 &&
        validation.specific_work > 0.0 &&
        validation.shaft_power > 0.0 &&
        validation.torque > 0.0 &&
        total_entropy >= 0.0;

    validation.operating_limits =
        OperatingLimitStatus::NotEvaluated;

    validation.message = independent
        ? "Coupled diffusion independently validated"
        : "Coupled diffusion independent validation failed";

    return {
        solved,
        validation,
        incidence,
        coupled.diffusion,
        Tt1,
        operating.entropy_generation,
        incidence_entropy,
        coupled.diffusion_entropy_generation,
        total_entropy,
        coupled.residual_norm,
        independent
    };
}

} // namespace turbocompgpu
