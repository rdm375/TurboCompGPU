#pragma once

#include <turbocompgpu/operating_point.hpp>
#include <turbocompgpu/vaneless_diffuser.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>

namespace turbocompgpu {

// Stage 8B.5e: reduced lossless static-pressure recovery boundary.
//
// Valid for r3 > r2, 0 < A3 < A2, and the positive-whirl
// branch of the current impeller velocity triangle.
//
// The diffuser equations are eliminated analytically.
// Each scalar residual evaluation solves only the impeller.

enum class RecoveryBoundaryStatus {
    Converged,
    InvalidInput,
    ImpellerFailed,
    NotBracketed,
    MaximumIterations
};

struct RecoveryBoundaryOptions {
    double minimum_mass_flow = 0.25;
    double maximum_mass_flow = 25.0;
    double scan_increment = 0.05;
    double relative_mass_tolerance = 1.0e-10;
    double relative_residual_tolerance = 1.0e-9;
    int maximum_iterations = 100;
};

struct RecoveryBoundaryResult {
    RecoveryBoundaryStatus status =
        RecoveryBoundaryStatus::InvalidInput;

    double mass_flow = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    double analytical_meridional_velocity =
        std::numeric_limits<double>::quiet_NaN();
    double outlet_density = std::numeric_limits<double>::quiet_NaN();
    double bracket_width = std::numeric_limits<double>::quiet_NaN();

    int impeller_evaluations = 0;
    int impeller_residual_evaluations = 0;
    int iterations = 0;

    std::string message;

    bool converged() const noexcept {
        return status == RecoveryBoundaryStatus::Converged;
    }
};

inline RecoveryBoundaryResult solve_recovery_boundary(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& impeller_geometry,
    const VanelessDiffuserGeometry& diffuser_geometry,
    OperatingPointRequest request,
    const NonlinearResidualScales& scales,
    const NonlinearSolverOptions& solver_options = {},
    const RecoveryBoundaryOptions& boundary_options = {},
    double validation_tolerance = 1.0e-8)
{
    RecoveryBoundaryResult result;

    const double r2 = impeller_geometry.exit.diameter / 2.0;
    const double r3 = diffuser_geometry.outlet_radius;

    const double area2 = effective_flow_area(
        impeller_exit_area(impeller_geometry.exit),
        impeller_geometry.exit.blockage);

    const double area3 =
        2.0 * std::numbers::pi *
        r3 * diffuser_geometry.outlet_width *
        (1.0 - diffuser_geometry.outlet_blockage);

    const double sigma = request.operating.slip_factor;
    const double rpm = request.operating.rotational_speed_rpm;

    const double beta =
        impeller_geometry.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const auto& settings = boundary_options;

    if (!(std::isfinite(r2) && std::isfinite(r3) &&
          std::isfinite(area2) && std::isfinite(area3) &&
          std::isfinite(sigma) && std::isfinite(rpm) &&
          std::isfinite(beta) &&
          r2 > 0.0 && r3 > r2 &&
          std::isfinite(diffuser_geometry.inlet_radius) &&
          std::abs(diffuser_geometry.inlet_radius - r2) <=
              1.0e-10 * r2 &&
          area3 > 0.0 && area3 < area2 &&
          diffuser_geometry.outlet_width > 0.0 &&
          diffuser_geometry.outlet_blockage >= 0.0 &&
          diffuser_geometry.outlet_blockage < 1.0 &&
          sigma > 0.0 && rpm > 0.0 &&
          beta > 0.0 && beta < std::numbers::pi / 2.0 &&
          std::isfinite(settings.minimum_mass_flow) &&
          std::isfinite(settings.maximum_mass_flow) &&
          std::isfinite(settings.scan_increment) &&
          settings.minimum_mass_flow > 0.0 &&
          settings.maximum_mass_flow > settings.minimum_mass_flow &&
          settings.scan_increment > 0.0 &&
          std::isfinite(settings.relative_mass_tolerance) &&
          settings.relative_mass_tolerance > 0.0 &&
          std::isfinite(settings.relative_residual_tolerance) &&
          settings.relative_residual_tolerance > 0.0 &&
          settings.maximum_iterations > 0)) {
        result.message = "Invalid analytical recovery configuration";
        return result;
    }

    const double alpha = area2 / area3;
    const double lambda = r3 / r2;

    const double K = std::sqrt(
        (alpha * alpha - 1.0) /
        (1.0 - 1.0 / (lambda * lambda)));

    const double omega =
        2.0 * std::numbers::pi * rpm / 60.0;

    const double vm = sigma * omega * r2 /
        (K + sigma / std::tan(beta));

    if (!(std::isfinite(vm) && vm > 0.0)) {
        result.message = "Invalid analytical meridional velocity";
        return result;
    }

    result.analytical_meridional_velocity = vm;

    struct Point {
        double mass = 0.0;
        double residual = 0.0;
        double density = 0.0;
        bool valid = false;
    };

    auto evaluate_point = [&](double mass) -> Point {
        request.operating.mass_flow = mass;

        const auto impeller = evaluate_operating_point(
            model, impeller_geometry, request, scales,
            solver_options, validation_tolerance);

        ++result.impeller_evaluations;
        result.impeller_residual_evaluations += impeller.evaluations;

        if (!(impeller.converged &&
              impeller.independently_validated &&
              impeller.basic_physics_valid &&
              std::isfinite(impeller.outlet_density) &&
              impeller.outlet_density > 0.0)) {
            return {};
        }

        const double residual =
            mass - impeller.outlet_density * area2 * vm;

        return {
            mass,
            residual,
            impeller.outlet_density,
            std::isfinite(residual)
        };
    };

    auto accept = [&](const Point& point, double width) {
        result.status = RecoveryBoundaryStatus::Converged;
        result.mass_flow = point.mass;
        result.residual = point.residual;
        result.outlet_density = point.density;
        result.bracket_width = width;
        result.message = "Reduced recovery boundary converged";
    };

    auto residual_small = [&](const Point& point) {
        return std::abs(point.residual) <=
            settings.relative_residual_tolerance *
            std::max(1.0, point.mass);
    };

    Point previous;
    bool have_previous = false;
    bool found_bracket = false;
    Point lower, upper;

    // Invalid intermediate evaluations break bracket continuity.
    const int scan_count = static_cast<int>(std::ceil(
        (settings.maximum_mass_flow -
         settings.minimum_mass_flow) / settings.scan_increment));

    for (int i = 0; i <= scan_count; ++i) {
        const double mass = std::min(
            settings.maximum_mass_flow,
            settings.minimum_mass_flow +
                i * settings.scan_increment);

        const Point current = evaluate_point(mass);

        if (!current.valid) {
            have_previous = false;
            continue;
        }

        if (residual_small(current)) {
            accept(current, 0.0);
            return result;
        }

        if (have_previous &&
            std::signbit(previous.residual) !=
            std::signbit(current.residual)) {
            lower = previous;
            upper = current;
            found_bracket = true;
            break;
        }

        previous = current;
        have_previous = true;
    }

    if (!found_bracket) {
        result.status = RecoveryBoundaryStatus::NotBracketed;
        result.message =
            "No continuous sign-changing bracket found";
        return result;
    }

    for (int iteration = 1;
         iteration <= settings.maximum_iterations;
         ++iteration) {

        result.iterations = iteration;

        const double midpoint =
            lower.mass + 0.5 * (upper.mass - lower.mass);

        if (midpoint == lower.mass || midpoint == upper.mass)
            break;

        const Point current = evaluate_point(midpoint);

        if (!current.valid) {
            result.status = RecoveryBoundaryStatus::ImpellerFailed;
            result.message =
                "Impeller evaluation failed inside root bracket";
            return result;
        }

        if (std::signbit(current.residual) ==
            std::signbit(lower.residual)) {
            lower = current;
        } else {
            upper = current;
        }

        const double width = upper.mass - lower.mass;

        if (width <= settings.relative_mass_tolerance *
                     std::max(1.0, midpoint)) {

            const Point final_point = evaluate_point(
                lower.mass + 0.5 * width);

            if (!final_point.valid) {
                result.status = RecoveryBoundaryStatus::ImpellerFailed;
                result.message = "Final impeller evaluation failed";
                return result;
            }

            if (residual_small(final_point)) {
                accept(final_point, width);
                return result;
            }

            result.mass_flow = final_point.mass;
            result.residual = final_point.residual;
            result.outlet_density = final_point.density;
            result.bracket_width = width;
            result.status = RecoveryBoundaryStatus::MaximumIterations;
            result.message =
                "Mass bracket converged but residual tolerance failed";
            return result;
        }
    }

    result.status = RecoveryBoundaryStatus::MaximumIterations;
    result.bracket_width = upper.mass - lower.mass;
    result.message = "Root refinement did not converge";
    return result;
}

} // namespace turbocompgpu
