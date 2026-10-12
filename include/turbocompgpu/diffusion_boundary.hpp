#pragma once

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace turbocompgpu {

struct DiffusionBoundaryOptions {
    double flow_tolerance = 1.0e-6;
    double diffusion_tolerance = 1.0e-8;
    int max_iterations = 80;
};

enum class DiffusionBoundaryStatus {
    converged,
    not_bracketed,
    maximum_iterations
};

struct DiffusionBoundaryResult {
    DiffusionBoundaryStatus status =
        DiffusionBoundaryStatus::not_bracketed;

    double critical_mass_flow =
        std::numeric_limits<double>::quiet_NaN();

    double diffusion_residual =
        std::numeric_limits<double>::quiet_NaN();

    double lower_mass_flow = 0.0;
    double upper_mass_flow = 0.0;

    int iterations = 0;
    int evaluations = 0;

    [[nodiscard]] bool converged() const noexcept {
        return status == DiffusionBoundaryStatus::converged;
    }
};

template<class Evaluator>
DiffusionBoundaryResult solve_diffusion_boundary(
    Evaluator&& evaluate_diffusion_ratio,
    double lower_mass_flow,
    double upper_mass_flow,
    double critical_diffusion_ratio,
    const DiffusionBoundaryOptions& options = {})
{
    if (!(std::isfinite(lower_mass_flow) &&
          std::isfinite(upper_mass_flow) &&
          lower_mass_flow > 0.0 &&
          upper_mass_flow > lower_mass_flow &&
          std::isfinite(critical_diffusion_ratio) &&
          critical_diffusion_ratio > 0.0 &&
          std::isfinite(options.flow_tolerance) &&
          options.flow_tolerance > 0.0 &&
          std::isfinite(options.diffusion_tolerance) &&
          options.diffusion_tolerance > 0.0 &&
          options.max_iterations > 0)) {
        throw std::invalid_argument(
            "Invalid diffusion boundary arguments");
    }

    DiffusionBoundaryResult result;
    result.lower_mass_flow = lower_mass_flow;
    result.upper_mass_flow = upper_mass_flow;

    auto residual = [&](double flow) {
        const double ratio = evaluate_diffusion_ratio(flow);
        ++result.evaluations;

        if (!std::isfinite(ratio)) {
            throw std::runtime_error(
                "Nonfinite diffusion ratio");
        }

        return ratio - critical_diffusion_ratio;
    };

    double a = lower_mass_flow;
    double b = upper_mass_flow;
    double fa = residual(a);
    double fb = residual(b);

    // An endpoint root is an exact zero-width bracket.
    if (fa == 0.0) {
        result.status = DiffusionBoundaryStatus::converged;
        result.critical_mass_flow = a;
        result.diffusion_residual = fa;
        result.upper_mass_flow = a;
        return result;
    }

    if (fb == 0.0) {
        result.status = DiffusionBoundaryStatus::converged;
        result.critical_mass_flow = b;
        result.diffusion_residual = fb;
        result.lower_mass_flow = b;
        return result;
    }

    if (std::signbit(fa) == std::signbit(fb)) {
        return result;
    }

    for (int iteration = 1;
         iteration <= options.max_iterations;
         ++iteration) {

        const double midpoint = std::midpoint(a, b);

        if (midpoint == a || midpoint == b) {
            break;
        }

        const double fm = residual(midpoint);
        result.iterations = iteration;

        if (fm == 0.0) {
            a = midpoint;
            b = midpoint;
            fa = fm;
            fb = fm;
        } else if (std::signbit(fm) == std::signbit(fa)) {
            a = midpoint;
            fa = fm;
        } else {
            b = midpoint;
            fb = fm;
        }

        result.lower_mass_flow = a;
        result.upper_mass_flow = b;

        if ((b - a) <= options.flow_tolerance) {
            const double root = std::midpoint(a, b);
            const double fr = residual(root);

            result.critical_mass_flow = root;
            result.diffusion_residual = fr;

            if (std::abs(fr) <= options.diffusion_tolerance) {
                result.status =
                    DiffusionBoundaryStatus::converged;
                return result;
            }
        }
    }

    result.status = DiffusionBoundaryStatus::maximum_iterations;
    result.lower_mass_flow = a;
    result.upper_mass_flow = b;
    result.critical_mass_flow = std::midpoint(a, b);
    result.diffusion_residual =
        residual(result.critical_mass_flow);

    return result;
}

} // namespace turbocompgpu
