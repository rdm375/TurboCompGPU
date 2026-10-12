#include <turbocompgpu/nonlinear_solver.hpp>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

struct SyntheticResidual {
    double normalized_energy;
    double normalized_entropy;
    double residual_norm;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

SyntheticResidual residual(double f1, double f2)
{
    return {f1, f2, std::hypot(f1, f2)};
}

using Evaluator =
    std::function<SyntheticResidual(const NonlinearUnknowns&)>;

using Result = BasicNonlinearSolverResult<SyntheticResidual>;

const NonlinearUnknowns initial {300.0, 0.0};

Result solve(
    const Evaluator& evaluator,
    const NonlinearSolverOptions& options = {},
    const NonlinearUnknowns& guess = initial)
{
    return solve_nonlinear_compressor_with_evaluator(
        evaluator, guess, options);
}

void check(
    const std::string& name,
    const Result& result,
    NonlinearSolveStatus expected,
    int expected_iterations,
    int expected_jacobians,
    int expected_rejections,
    bool expect_solution)
{
    require(
        result.status == expected,
        name + ": unexpected status " +
            nonlinear_status_name(result.status));

    require(
        result.converged() ==
            (expected == NonlinearSolveStatus::Converged),
        name + ": inconsistent converged()");

    require(
        result.solution.has_value() == expect_solution,
        name + ": unexpected solution presence");

    require(
        result.iterations == expected_iterations,
        name + ": unexpected iteration count");

    require(
        result.jacobian_evaluations == expected_jacobians,
        name + ": unexpected Jacobian count");

    require(
        result.rejected_steps == expected_rejections,
        name + ": unexpected rejected-step count");

    require(
        !result.message.empty(),
        name + ": missing diagnostic message");

    std::cout
        << name
        << " status=" << nonlinear_status_name(result.status)
        << " iterations=" << result.iterations
        << " evaluations=" << result.residual_evaluations
        << " jacobians=" << result.jacobian_evaluations
        << " rejected=" << result.rejected_steps
        << '\n';
}

void test_invalid_input()
{
    int calls = 0;

    const Evaluator evaluator =
        [&](const NonlinearUnknowns&) {
            ++calls;
            return residual(1.0, 1.0);
        };

    const NonlinearUnknowns invalid {
        -300.0, 0.0
    };

    const auto result = solve(evaluator, {}, invalid);

    check(
        "InvalidInput",
        result,
        NonlinearSolveStatus::InvalidInput,
        0, 0, 0, false);

    require(calls == 0, "Invalid input evaluated residual");
    require(
        result.residual_evaluations == 0,
        "Invalid input incremented evaluation count");

    NonlinearSolverOptions options;
    options.difference_step = 0.0;

    const auto bad_options = solve(evaluator, options);

    check(
        "InvalidOptions",
        bad_options,
        NonlinearSolveStatus::InvalidInput,
        0, 0, 0, false);

    require(calls == 0, "Invalid options evaluated residual");
}

void test_initial_evaluation_failed()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns&) -> SyntheticResidual {
            throw std::runtime_error("Synthetic EOS failure");
        };

    const auto result = solve(evaluator);

    check(
        "InitialEvaluationFailed",
        result,
        NonlinearSolveStatus::InitialEvaluationFailed,
        0, 0, 0, false);

    require(
        result.residual_evaluations == 1,
        "Initial failure evaluation count incorrect");
}

void test_jacobian_failed()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns& x) {
            // Only the initial state is admissible.
            // Both perturbations fail for the first column.
            if (x.outlet_temperature != 300.0 ||
                x.log_pressure != 0.0) {
                throw std::runtime_error(
                    "Synthetic EOS domain boundary");
            }

            return residual(1.0, 1.0);
        };

    const auto result = solve(evaluator);

    check(
        "JacobianFailed",
        result,
        NonlinearSolveStatus::JacobianFailed,
        0, 1, 0, true);

    require(
        result.residual_evaluations == 3,
        "Jacobian failure evaluation count incorrect");
}

void test_singular_jacobian()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns& x) {
            const double y =
                x.outlet_temperature / 300.0 - 1.0;

            const double z = x.log_pressure;

            // Two identical residual components:
            // the Jacobian has rank one.
            const double f = y + z + 1.0;

            return residual(f, f);
        };

    const auto result = solve(evaluator);

    check(
        "SingularJacobian",
        result,
        NonlinearSolveStatus::SingularJacobian,
        0, 1, 0, true);

    require(
        result.residual_evaluations == 5,
        "Singular Jacobian evaluation count incorrect");
}

void test_line_search_failed()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns& x) {
            const double y =
                x.outlet_temperature / 300.0 - 1.0;

            const double z = x.log_pressure;

            // Locally well-conditioned at the initial point.
            // Newton's proposed full step is rejected because
            // every noninitial trial state has a larger norm.
            //
            // The discontinuity is deliberate: this is a
            // controlled line-search failure, not a physical EOS.
            if (y == 0.0 && z == 0.0)
                return residual(1.0, 1.0);

            if (std::abs(y) < 1.0e-3 &&
                std::abs(z) < 1.0e-3)
                return residual(1.0 + y, 1.0 + z);

            return residual(10.0, 10.0);
        };

    NonlinearSolverOptions options;
    options.maximum_backtracks = 0;

    const auto result = solve(evaluator, options);

    check(
        "LineSearchFailed",
        result,
        NonlinearSolveStatus::LineSearchFailed,
        0, 1, 1, true);

    require(
        result.residual_evaluations == 6,
        "Line search failure evaluation count incorrect");
}

void test_maximum_iterations()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns& x) {
            return residual(
                x.outlet_temperature / 300.0,
                x.log_pressure + 1.0);
        };

    NonlinearSolverOptions options;
    options.maximum_iterations = 0;

    const auto result = solve(evaluator, options);

    check(
        "MaximumIterations",
        result,
        NonlinearSolveStatus::MaximumIterations,
        0, 0, 0, true);

    require(
        result.residual_evaluations == 1,
        "Maximum iterations evaluation count incorrect");
}

void test_converged()
{
    const Evaluator evaluator =
        [](const NonlinearUnknowns& x) {
            return residual(
                x.outlet_temperature / 300.0 - 1.0,
                x.log_pressure);
        };

    const auto result = solve(evaluator);

    check(
        "Converged",
        result,
        NonlinearSolveStatus::Converged,
        0, 0, 0, true);

    require(
        result.final_residual_norm == 0.0,
        "Exact solution has nonzero residual");

    require(
        result.residual_evaluations == 1,
        "Converged evaluation count incorrect");
}

} // namespace

int main()
{
    try {
        test_invalid_input();
        test_initial_evaluation_failed();
        test_jacobian_failed();
        test_singular_jacobian();
        test_line_search_failed();
        test_maximum_iterations();
        test_converged();

        std::cout
            << "PASS: All nonlinear solver failure paths\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
