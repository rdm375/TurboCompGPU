#pragma once

#include <turbocompgpu/nonlinear_residual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace turbocompgpu {

enum class NonlinearSolveStatus {
    Converged,
    InvalidInput,
    InitialEvaluationFailed,
    JacobianFailed,
    SingularJacobian,
    LineSearchFailed,
    MaximumIterations
};

inline const char* nonlinear_status_name(
    NonlinearSolveStatus status) noexcept
{
    switch (status) {
    case NonlinearSolveStatus::Converged:
        return "Converged";
    case NonlinearSolveStatus::InvalidInput:
        return "InvalidInput";
    case NonlinearSolveStatus::InitialEvaluationFailed:
        return "InitialEvaluationFailed";
    case NonlinearSolveStatus::JacobianFailed:
        return "JacobianFailed";
    case NonlinearSolveStatus::SingularJacobian:
        return "SingularJacobian";
    case NonlinearSolveStatus::LineSearchFailed:
        return "LineSearchFailed";
    case NonlinearSolveStatus::MaximumIterations:
        return "MaximumIterations";
    }
    return "Unknown";
}

struct NonlinearSolverOptions {
    int maximum_iterations = 40;
    int maximum_backtracks = 24;

    double residual_tolerance = 1.0e-10;
    double difference_step = 1.0e-5;

    double temperature_scale = 300.0;

    double minimum_step = 1.0e-8;
    double jacobian_singularity_tolerance = 1.0e-12;
};

struct NonlinearSolverResult {
    NonlinearSolveStatus status =
        NonlinearSolveStatus::InvalidInput;

    std::optional<NonlinearResidualResult> solution;

    int iterations = 0;
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
    int rejected_steps = 0;

    double final_residual_norm =
        std::numeric_limits<double>::infinity();

    std::string message;

    bool converged() const noexcept
    {
        return status == NonlinearSolveStatus::Converged;
    }
};

namespace nonlinear_detail {

using Vector = std::array<double, 2>;

inline Vector residual_vector(
    const NonlinearResidualResult& r)
{
    return {r.normalized_energy, r.normalized_entropy};
}

inline bool valid_options(const NonlinearSolverOptions& o)
{
    return o.maximum_iterations >= 0 &&
           o.maximum_backtracks >= 0 &&
           std::isfinite(o.residual_tolerance) &&
           o.residual_tolerance > 0.0 &&
           std::isfinite(o.difference_step) &&
           o.difference_step > 0.0 &&
           std::isfinite(o.temperature_scale) &&
           o.temperature_scale > 0.0 &&
           std::isfinite(o.minimum_step) &&
           o.minimum_step > 0.0 &&
           o.minimum_step <= 1.0 &&
           std::isfinite(o.jacobian_singularity_tolerance) &&
           o.jacobian_singularity_tolerance > 0.0;
}

inline bool solve_2x2(
    const std::array<Vector, 2>& columns,
    const Vector& rhs,
    double singularity_tolerance,
    Vector& step)
{
    // Matrix stored as two columns.
    const double a = columns[0][0];
    const double c = columns[0][1];
    const double b = columns[1][0];
    const double d = columns[1][1];

    const double norm0 = std::hypot(a, c);
    const double norm1 = std::hypot(b, d);

    if (!(norm0 > 0.0 && norm1 > 0.0) ||
        !std::isfinite(norm0) ||
        !std::isfinite(norm1)) {
        return false;
    }

    // Column scaling avoids unnecessary overflow in
    // the determinant and makes the singularity check
    // independent of column magnitudes.
    const double as = a / norm0;
    const double cs = c / norm0;
    const double bs = b / norm1;
    const double ds = d / norm1;

    const double det = as * ds - bs * cs;

    if (!std::isfinite(det) ||
        std::abs(det) <= singularity_tolerance) {
        return false;
    }

    const double rhs0 = rhs[0];
    const double rhs1 = rhs[1];

    step[0] =
        (ds * rhs0 - bs * rhs1) / det / norm0;

    step[1] =
        (-cs * rhs0 + as * rhs1) / det / norm1;

    return std::isfinite(step[0]) &&
           std::isfinite(step[1]);
}

} // namespace nonlinear_detail

inline NonlinearSolverResult solve_nonlinear_compressor_with_evaluator(
    const std::function<NonlinearResidualResult(
        const NonlinearUnknowns&)>& evaluator,
    const NonlinearUnknowns& initial_guess,
    const NonlinearSolverOptions& options = {})
{
    using nonlinear_detail::Vector;

    NonlinearSolverResult result;

    if (!nonlinear_detail::valid_options(options) ||
        !std::isfinite(initial_guess.outlet_temperature) ||
        initial_guess.outlet_temperature <= 0.0 ||
        !std::isfinite(initial_guess.log_pressure)) {
        result.message = "Invalid solver options or initial guess";
        return result;
    }

    // Dimensionless solver coordinates.
    Vector y {
        initial_guess.outlet_temperature /
            options.temperature_scale,
        initial_guess.log_pressure
    };

    auto evaluate = [&](const Vector& point)
        -> std::optional<NonlinearResidualResult>
    {
        ++result.residual_evaluations;

        if (!std::isfinite(point[0]) ||
            !std::isfinite(point[1]) ||
            point[0] <= 0.0) {
            return std::nullopt;
        }

        try {
            auto r = evaluator({
                point[0] * options.temperature_scale,
                point[1]
            });

            if (!std::isfinite(r.residual_norm)) {
                return std::nullopt;
            }

            return r;
        }
        catch (const std::exception&) {
            // Invalid thermodynamic trial states are
            // rejected during differencing or line search.
            return std::nullopt;
        }
    };

    auto current = evaluate(y);

    if (!current) {
        result.status =
            NonlinearSolveStatus::InitialEvaluationFailed;
        result.message =
            "Initial state could not be evaluated";
        return result;
    }

    auto record_current = [&]() {
        result.solution = *current;
        result.final_residual_norm =
            current->residual_norm;
    };

    record_current();

    for (int iteration = 0;
         iteration <= options.maximum_iterations;
         ++iteration) {

        if (current->residual_norm <=
            options.residual_tolerance) {
            result.status = NonlinearSolveStatus::Converged;
            result.message =
                "Normalized residual tolerance satisfied";
            record_current();
            return result;
        }

        if (iteration == options.maximum_iterations) {
            break;
        }

        ++result.jacobian_evaluations;

        const Vector f0 =
            nonlinear_detail::residual_vector(*current);

        std::array<Vector, 2> columns {};

        for (int j = 0; j < 2; ++j) {
            const double h =
                options.difference_step *
                std::max(1.0, std::abs(y[j]));

            Vector yp = y;
            Vector ym = y;

            yp[j] += h;
            ym[j] -= h;

            const auto plus = evaluate(yp);
            const auto minus = evaluate(ym);

            Vector derivative {};

            if (plus && minus) {
                const auto fp =
                    nonlinear_detail::residual_vector(*plus);
                const auto fm =
                    nonlinear_detail::residual_vector(*minus);

                const double denominator =
                    yp[j] - ym[j];

                for (int i = 0; i < 2; ++i) {
                    derivative[i] =
                        (fp[i] - fm[i]) / denominator;
                }
            }
            else if (plus) {
                const auto fp =
                    nonlinear_detail::residual_vector(*plus);

                const double denominator =
                    yp[j] - y[j];

                for (int i = 0; i < 2; ++i) {
                    derivative[i] =
                        (fp[i] - f0[i]) / denominator;
                }
            }
            else if (minus) {
                const auto fm =
                    nonlinear_detail::residual_vector(*minus);

                const double denominator =
                    y[j] - ym[j];

                for (int i = 0; i < 2; ++i) {
                    derivative[i] =
                        (f0[i] - fm[i]) / denominator;
                }
            }
            else {
                result.status =
                    NonlinearSolveStatus::JacobianFailed;
                result.message =
                    "No valid finite-difference perturbation";
                record_current();
                return result;
            }

            if (!std::isfinite(derivative[0]) ||
                !std::isfinite(derivative[1])) {
                result.status =
                    NonlinearSolveStatus::JacobianFailed;
                result.message =
                    "Non-finite Jacobian derivative";
                record_current();
                return result;
            }

            columns[j] = derivative;
        }

        Vector step {};

        if (!nonlinear_detail::solve_2x2(
                columns,
                {-f0[0], -f0[1]},
                options.jacobian_singularity_tolerance,
                step)) {
            result.status =
                NonlinearSolveStatus::SingularJacobian;
            result.message =
                "Singular or nearly singular Jacobian";
            record_current();
            return result;
        }

        // Armijo backtracking on the Euclidean norm
        // of the dimensionless residual vector.
        //
        // For an exact Newton direction:
        //
        //   J * step = -F
        //
        // so the directional derivative of ||F|| is
        // -||F|| when F != 0.

        constexpr double armijo = 1.0e-4;

        double alpha = 1.0;
        bool accepted = false;

        for (int backtrack = 0;
             backtrack <= options.maximum_backtracks;
             ++backtrack) {

            if (alpha < options.minimum_step) {
                break;
            }

            const Vector candidate_y {
                y[0] + alpha * step[0],
                y[1] + alpha * step[1]
            };

            auto candidate = evaluate(candidate_y);

            if (candidate &&
                candidate->residual_norm <=
                    (1.0 - armijo * alpha) *
                    current->residual_norm) {

                y = candidate_y;
                current = std::move(candidate);
                accepted = true;

                ++result.iterations;
                record_current();
                break;
            }

            ++result.rejected_steps;
            alpha *= 0.5;
        }

        if (!accepted) {
            result.status =
                NonlinearSolveStatus::LineSearchFailed;
            result.message =
                "No admissible residual-decreasing Newton step";
            record_current();
            return result;
        }
    }

    result.status =
        NonlinearSolveStatus::MaximumIterations;
    result.message =
        "Maximum Newton iterations reached";

    record_current();
    return result;
}

inline NonlinearSolverResult solve_nonlinear_compressor(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const NonlinearOperatingPoint& operating,
    const NonlinearUnknowns& initial_guess,
    const NonlinearResidualScales& scales,
    const NonlinearSolverOptions& options = {})
{
    const auto evaluator =
        [&](const NonlinearUnknowns& unknowns) {
            return evaluate_nonlinear_residual(
                model,
                geometry,
                operating,
                unknowns,
                scales);
        };

    return solve_nonlinear_compressor_with_evaluator(
        evaluator,
        initial_guess,
        options);
}

} // namespace turbocompgpu
