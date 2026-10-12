#pragma once

#include <turbocompgpu/choking_continuation.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace turbocompgpu {

// Stage 6H.1: matched physical-geometry and constant-area
// blade-angle choking comparison.
//
// The constant-area geometry changes passage height to
// preserve the reference effective throat area.
//
// This is a controlled geometric comparison, not an
// additive decomposition of aerodynamic losses.

struct ChokingAreaComparisonPoint {
    double blade_angle_deg;

    double physical_height;
    double constant_area_height;

    double physical_area;
    double constant_area;

    double physical_mass_flow;
    double constant_area_mass_flow;

    double physical_mass_flux;
    double constant_area_mass_flux;

    double physical_tangential_velocity;
    double constant_area_tangential_velocity;

    double physical_mass_residual;
    double constant_area_mass_residual;

    double physical_rothalpy_residual;
    double constant_area_rothalpy_residual;

    double mass_flow_difference;
    double mass_flow_ratio;
};

inline ThroatGeometry
make_constant_area_throat(
    const ThroatGeometry& reference,
    double blade_angle_deg)
{
    const auto reference_result =
        evaluate_throat_geometry(reference);

    ThroatGeometry adjusted = reference;
    adjusted.blade_angle_deg = blade_angle_deg;

    // Validate the new angle and passage opening using
    // the original height before calculating its replacement.
    const auto provisional =
        evaluate_throat_geometry(adjusted);

    const double area_per_height =
        static_cast<double>(adjusted.blade_count) *
        provisional.open_width *
        (1.0 - adjusted.blockage);

    if (!(std::isfinite(area_per_height) &&
          area_per_height > 0.0)) {
        throw std::runtime_error(
            "Invalid constant-area throat denominator");
    }

    adjusted.passage_height =
        reference_result.effective_area / area_per_height;

    const auto verified =
        evaluate_throat_geometry(adjusted);

    const double tolerance =
        1e-12 * std::max(
            1.0, reference_result.effective_area);

    if (std::abs(
            verified.effective_area -
            reference_result.effective_area) > tolerance) {
        throw std::runtime_error(
            "Constant-area throat verification failed");
    }

    return adjusted;
}

template<class Model>
std::vector<ChokingAreaComparisonPoint>
evaluate_choking_area_comparison(
    const Model& model,
    const ImpellerSystemGeometry& geometry,
    const PassageTurningChokingInputs& baseline,
    const std::vector<double>& blade_angles)
{
    std::vector<ChokingAreaComparisonPoint> points;
    points.reserve(blade_angles.size());

    for (double angle : blade_angles) {
        auto physical_input = baseline;
        physical_input.throat_geometry.blade_angle_deg =
            angle;

        auto constant_input = baseline;
        constant_input.throat_geometry =
            make_constant_area_throat(
                baseline.throat_geometry, angle);

        const auto physical_geometry =
            evaluate_throat_geometry(
                physical_input.throat_geometry);

        const auto constant_geometry =
            evaluate_throat_geometry(
                constant_input.throat_geometry);

        const auto physical =
            evaluate_passage_turning_choking(
                model, model, geometry, physical_input);

        const auto constant =
            evaluate_passage_turning_choking(
                model, model, geometry, constant_input);

        const double mp =
            physical.choking.coupled.critical_mass_flow;

        const double mc =
            constant.choking.coupled.critical_mass_flow;

        if (!(std::isfinite(mp) && mp > 0.0 &&
              std::isfinite(mc) && mc > 0.0)) {
            throw std::runtime_error(
                "Invalid choking mass flow in area comparison");
        }

        points.push_back({
            angle,
            physical_input.throat_geometry.passage_height,
            constant_input.throat_geometry.passage_height,
            physical_geometry.effective_area,
            constant_geometry.effective_area,
            mp,
            mc,
            mp / physical_geometry.effective_area,
            mc / constant_geometry.effective_area,
            physical.predicted_tangential_velocity,
            constant.predicted_tangential_velocity,
            physical.mass_flow_residual,
            constant.mass_flow_residual,
            physical.choking.rothalpy_residual,
            constant.choking.rothalpy_residual,
            mp - mc,
            mp / mc
        });
    }

    return points;
}

} // namespace turbocompgpu
