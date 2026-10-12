#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/nonlinear_solver.hpp>

#include <concepts>
#include <optional>
#include <type_traits>
#include <utility>

using namespace turbocompgpu;

template <typename T>
concept HasFullSystem = requires(T value) {
    value.full_system;
};

static_assert(
    HasFullSystem<NonlinearResidualResult>);

static_assert(
    !HasFullSystem<FastNonlinearResidualResult>);

static_assert(
    !std::is_same_v<
        NonlinearResidualResult,
        FastNonlinearResidualResult>);

static_assert(
    std::is_same_v<
        decltype(std::declval<NonlinearSolverResult>().solution),
        std::optional<NonlinearResidualResult>>);

static_assert(
    std::is_same_v<
        decltype(std::declval<FastNonlinearSolverResult>().solution),
        std::optional<FastNonlinearResidualResult>>);

static_assert(
    std::is_same_v<
        decltype(solve_fast_nonlinear_compressor(
            std::declval<const FastResidualContext&>(),
            std::declval<const NonlinearUnknowns&>())),
        FastNonlinearSolverResult>);

static_assert(
    std::is_same_v<
        decltype(solve_nonlinear_compressor(
            std::declval<const thermomodel::Model&>(),
            std::declval<const ImpellerSystemGeometry&>(),
            std::declval<const NonlinearOperatingPoint&>(),
            std::declval<const NonlinearUnknowns&>(),
            std::declval<const NonlinearResidualScales&>())),
        NonlinearSolverResult>);

int main()
{
    return 0;
}
