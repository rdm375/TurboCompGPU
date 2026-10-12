#pragma once

#include <turbocompgpu/fast_nonlinear_residual.hpp>
#include <turbocompgpu/nonlinear_solver.hpp>

namespace turbocompgpu {

using FastNonlinearSolverResult =
    BasicNonlinearSolverResult<FastNonlinearResidualResult>;

inline FastNonlinearSolverResult solve_fast_nonlinear_compressor(
    const FastResidualContext& context,
    const NonlinearUnknowns& initial_guess,
    const NonlinearSolverOptions& options = {})
{
    const auto evaluator =
        [&](const NonlinearUnknowns& unknowns) {
            return evaluate_fast_nonlinear_residual(
                context, unknowns);
        };

    return solve_nonlinear_compressor_with_evaluator(
        evaluator,
        initial_guess,
        options);
}

} // namespace turbocompgpu
