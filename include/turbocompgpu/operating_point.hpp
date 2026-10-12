#pragma once

#include <turbocompgpu/fast_nonlinear_solver.hpp>

#include <cmath>
#include <limits>
#include <string>

namespace turbocompgpu {

struct OperatingPointRequest {
    NonlinearOperatingPoint operating;
    NonlinearUnknowns initial_guess;
};

enum class OperatingLimitStatus {
    NotEvaluated,
    Admissible,
    Inadmissible
};

struct OperatingPointResult {
    NonlinearOperatingPoint operating;
    NonlinearSolveStatus status = NonlinearSolveStatus::InvalidInput;

    bool converged = false;
    bool independently_validated = false;
    bool basic_physics_valid = false;

    // Surge and choke limits are not modeled in M1.5.
    OperatingLimitStatus operating_limits =
        OperatingLimitStatus::NotEvaluated;

    double inlet_density =
        std::numeric_limits<double>::quiet_NaN();
    double outlet_density =
        std::numeric_limits<double>::quiet_NaN();

    double outlet_temperature = std::numeric_limits<double>::quiet_NaN();
    double outlet_pressure = std::numeric_limits<double>::quiet_NaN();
    double pressure_ratio = std::numeric_limits<double>::quiet_NaN();
    double specific_work = std::numeric_limits<double>::quiet_NaN();
    double shaft_power = std::numeric_limits<double>::quiet_NaN();
    double torque = std::numeric_limits<double>::quiet_NaN();
    double entropy_generation = std::numeric_limits<double>::quiet_NaN();
    double residual_norm = std::numeric_limits<double>::infinity();
    double validation_residual_norm = std::numeric_limits<double>::infinity();

    int iterations = 0;
    int evaluations = 0;
    int jacobians = 0;
    int rejected_steps = 0;

    std::string message;
};

// Independent validation here means agreement with the full
// residual evaluator. It does not establish surge/choke admissibility.
inline OperatingPointResult evaluate_operating_point(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const NonlinearSolverOptions& options = {},
    double validation_tolerance = 1.0e-8)
{
    OperatingPointResult result;
    result.operating = request.operating;

    if (!std::isfinite(validation_tolerance) ||
        validation_tolerance <= 0.0) {
        result.message = "Invalid validation tolerance";
        return result;
    }

    try {
        const auto context = make_fast_residual_context(
            model, geometry, request.operating, scales);

        const auto solved = solve_fast_nonlinear_compressor(
            context, request.initial_guess, options);

        result.status = solved.status;
        result.converged = solved.converged();
        result.iterations = solved.iterations;
        result.evaluations = solved.residual_evaluations;
        result.jacobians = solved.jacobian_evaluations;
        result.rejected_steps = solved.rejected_steps;
        result.residual_norm = solved.final_residual_norm;
        result.message = solved.message;

        if (!solved.solution)
            return result;

        const auto& solution = *solved.solution;

        result.outlet_temperature =
            solution.unknowns.outlet_temperature;
        result.outlet_pressure = solution.outlet_pressure;
        result.pressure_ratio =
            solution.outlet_pressure / request.operating.inlet_pressure;
        result.specific_work = solution.euler_specific_work;
        result.shaft_power =
            request.operating.mass_flow * solution.euler_specific_work;
        result.torque = solution.shaft_torque;
        result.entropy_generation = request.operating.entropy_generation;

        if (!result.converged)
            return result;

        // Independently reconstruct the complete physical system.
        const auto full = evaluate_nonlinear_residual(
            model, geometry, request.operating,
            solution.unknowns, scales);

        result.validation_residual_norm = full.residual_norm;

        result.independently_validated =
            std::isfinite(full.residual_norm) &&
            full.residual_norm <= validation_tolerance &&
            std::abs(full.euler_specific_work -
                     solution.euler_specific_work) <=
                validation_tolerance *
                std::max(1.0, std::abs(full.euler_specific_work)) &&
            std::abs(full.shaft_torque - solution.shaft_torque) <=
                validation_tolerance *
                std::max(1.0, std::abs(full.shaft_torque));

        if (!result.independently_validated) {
            result.message += "; independent validation failed";
            return result;
        }

        const auto outlet_state = model.evaluate_Tp(
            result.outlet_temperature,
            result.outlet_pressure);

        if (!outlet_state) {
            result.message += "; outlet state validation failed";
            return result;
        }

        result.inlet_density = context.inlet_state.rho;
        result.outlet_density = outlet_state.value.rho;

        const double omega = context.angular_speed;
        const double torque_power = result.torque * omega;

        const double power_scale =
            std::max(1.0, std::abs(result.shaft_power));

        const bool power_consistent =
            std::isfinite(torque_power) &&
            std::abs(result.shaft_power - torque_power)
                <= validation_tolerance * power_scale;

        result.basic_physics_valid =
            std::isfinite(request.operating.inlet_temperature) &&
            request.operating.inlet_temperature > 0.0 &&
            std::isfinite(request.operating.inlet_pressure) &&
            request.operating.inlet_pressure > 0.0 &&
            std::isfinite(result.outlet_temperature) &&
            result.outlet_temperature > 0.0 &&
            std::isfinite(result.outlet_pressure) &&
            result.outlet_pressure > 0.0 &&
            std::isfinite(result.inlet_density) &&
            result.inlet_density > 0.0 &&
            std::isfinite(result.outlet_density) &&
            result.outlet_density > 0.0 &&
            std::isfinite(result.pressure_ratio) &&
            result.pressure_ratio > 1.0 &&
            std::isfinite(result.specific_work) &&
            result.specific_work > 0.0 &&
            std::isfinite(result.shaft_power) &&
            result.shaft_power > 0.0 &&
            std::isfinite(result.torque) &&
            result.torque > 0.0 &&
            std::isfinite(result.entropy_generation) &&
            result.entropy_generation >= 0.0 &&
            power_consistent;

        if (!result.basic_physics_valid)
            result.message += "; basic physical checks failed";

        // Do not infer surge/choke admissibility from these checks.
        result.operating_limits =
            OperatingLimitStatus::NotEvaluated;
    }
    catch (const std::exception& error) {
        result.message +=
            std::string("; evaluation error: ") + error.what();
    }

    return result;
}

} // namespace turbocompgpu
