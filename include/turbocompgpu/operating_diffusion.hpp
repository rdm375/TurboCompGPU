#pragma once

#include <turbocompgpu/diffusion_loss.hpp>
#include <turbocompgpu/impeller_system.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct OperatingDiffusionConditions {
    double inlet_meridional_velocity;
    double outlet_meridional_velocity;

    double inlet_whirl_velocity;
    double outlet_whirl_velocity;

    double rotational_speed_rpm;
};

struct OperatingDiffusionResult {
    double inlet_blade_speed;
    double outlet_blade_speed;

    double inlet_relative_tangential_velocity;
    double outlet_relative_tangential_velocity;

    DiffusionLossResult loss;
};

inline OperatingDiffusionResult
evaluate_operating_diffusion_loss(
    const ImpellerSystemGeometry& geometry,
    const OperatingDiffusionConditions& operating,
    const DiffusionLossParameters& parameters)
{
    const double D1 =
        geometry.impeller.inlet_diameter;

    const double D2 =
        geometry.impeller.exit_diameter;

    if (!std::isfinite(D1) ||
        !std::isfinite(D2) ||
        D1 <= 0.0 ||
        D2 <= D1 ||
        !std::isfinite(operating.rotational_speed_rpm) ||
        operating.rotational_speed_rpm <= 0.0 ||
        !std::isfinite(operating.inlet_meridional_velocity) ||
        operating.inlet_meridional_velocity < 0.0 ||
        !std::isfinite(operating.outlet_meridional_velocity) ||
        operating.outlet_meridional_velocity < 0.0 ||
        !std::isfinite(operating.inlet_whirl_velocity) ||
        !std::isfinite(operating.outlet_whirl_velocity)) {
        throw std::invalid_argument(
            "Invalid operating-diffusion inputs"
        );
    }

    const double omega =
        2.0 * std::numbers::pi *
        operating.rotational_speed_rpm / 60.0;

    const double U1 = 0.5 * D1 * omega;
    const double U2 = 0.5 * D2 * omega;

    const double Wtheta1 =
        operating.inlet_whirl_velocity - U1;

    const double Wtheta2 =
        operating.outlet_whirl_velocity - U2;

    const double W1 = std::hypot(
        operating.inlet_meridional_velocity,
        Wtheta1
    );

    const double W2 = std::hypot(
        operating.outlet_meridional_velocity,
        Wtheta2
    );

    const auto loss = evaluate_diffusion_loss(
        W1, W2, parameters
    );

    return {
        U1,
        U2,
        Wtheta1,
        Wtheta2,
        loss
    };
}

} // namespace turbocompgpu
