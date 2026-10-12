#include <turbocompgpu/diffusion_boundary.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace turbocompgpu;

int main()
{
    // 1. Linear function with known root at m = 4.
    {
        const auto result = solve_diffusion_boundary(
            [](double m) { return 2.5 - 0.25 * m; },
            2.0, 6.0, 1.5);

        assert(result.converged());
        assert(std::abs(result.critical_mass_flow - 4.0) < 1e-6);
        assert(std::abs(result.diffusion_residual) < 1e-8);
        assert(result.upper_mass_flow -
               result.lower_mass_flow <= 1e-6);
    }

    // 2. Nonlinear function with known root at m = 3.
    {
        const auto result = solve_diffusion_boundary(
            [](double m) { return 1.5 + (9.0 - m*m) / 100.0; },
            2.0, 5.0, 1.5);

        assert(result.converged());
        assert(std::abs(result.critical_mass_flow - 3.0) < 1e-6);
    }

    // 3. No sign-changing bracket.
    {
        const auto result = solve_diffusion_boundary(
            [](double m) { return 2.0 + m / 100.0; },
            2.0, 8.0, 1.5);

        assert(!result.converged());
        assert(result.status ==
               DiffusionBoundaryStatus::not_bracketed);
        assert(result.evaluations == 2);
    }

    // 4. Exact endpoint root.
    {
        const auto result = solve_diffusion_boundary(
            [](double m) { return 2.5 - 0.25 * m; },
            4.0, 8.0, 1.5);

        assert(result.converged());
        assert(result.critical_mass_flow == 4.0);
        assert(result.diffusion_residual == 0.0);
    }

    // 5. Iteration limit must not produce false convergence.
    {
        DiffusionBoundaryOptions options;
        options.max_iterations = 1;

        const auto result = solve_diffusion_boundary(
            [](double m) { return 2.5 - 0.25 * m; },
            2.0, 8.0, 1.5, options);

        assert(!result.converged());
        assert(result.status ==
               DiffusionBoundaryStatus::maximum_iterations);
    }

    // 6. Invalid arguments must be rejected.
    {
        bool threw = false;

        try {
            (void)solve_diffusion_boundary(
                [](double m) { return m; },
                8.0, 2.0, 1.5);
        } catch (const std::invalid_argument&) {
            threw = true;
        }

        assert(threw);
    }

    // 7. Nonfinite evaluator results must be rejected.
    {
        bool threw = false;

        try {
            (void)solve_diffusion_boundary(
                [](double) {
                    return std::numeric_limits<double>::quiet_NaN();
                },
                2.0, 8.0, 1.5);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        assert(threw);
    }

    std::cout << "PASS: diffusion boundary root finder\n";
}
