#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct IncidenceLossParameters {
    double blade_inlet_angle_deg;
    double loss_coefficient;
};

struct IncidenceLossResult {
    double incidence_angle_deg;
    double relative_velocity;
    double specific_loss;
};

inline IncidenceLossResult evaluate_incidence_loss(
    double relative_flow_angle_deg,
    double relative_velocity,
    const IncidenceLossParameters& params)
{
    if (!std::isfinite(relative_flow_angle_deg) ||
        !std::isfinite(relative_velocity) ||
        relative_velocity < 0.0 ||
        !std::isfinite(params.blade_inlet_angle_deg) ||
        !std::isfinite(params.loss_coefficient) ||
        params.loss_coefficient < 0.0) {
        throw std::invalid_argument(
            "Invalid incidence-loss inputs"
        );
    }

    const double incidence_deg =
        relative_flow_angle_deg -
        params.blade_inlet_angle_deg;

    const double incidence_rad =
        incidence_deg *
        std::numbers::pi / 180.0;

    const double transverse_velocity =
        relative_velocity *
        std::sin(incidence_rad);

    const double specific_loss =
        0.5 *
        params.loss_coefficient *
        transverse_velocity *
        transverse_velocity;

    if (!std::isfinite(specific_loss)) {
        throw std::runtime_error(
            "Non-finite incidence loss"
        );
    }

    return {
        incidence_deg,
        relative_velocity,
        specific_loss
    };
}

} // namespace turbocompgpu
