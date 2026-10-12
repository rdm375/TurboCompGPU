#pragma once

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct DiffusionLossParameters {
    double loss_coefficient;
    double critical_diffusion_ratio;
};

struct DiffusionLossResult {
    double inlet_relative_velocity;
    double outlet_relative_velocity;
    double diffusion_ratio;
    double excess_diffusion;
    double specific_loss;
};

inline DiffusionLossResult evaluate_diffusion_loss(
    double inlet_relative_velocity,
    double outlet_relative_velocity,
    const DiffusionLossParameters& params)
{
    if (!std::isfinite(inlet_relative_velocity) ||
        inlet_relative_velocity < 0.0 ||
        !std::isfinite(outlet_relative_velocity) ||
        outlet_relative_velocity <= 0.0 ||
        !std::isfinite(params.loss_coefficient) ||
        params.loss_coefficient < 0.0 ||
        !std::isfinite(params.critical_diffusion_ratio) ||
        params.critical_diffusion_ratio < 1.0) {
        throw std::invalid_argument(
            "Invalid diffusion-loss inputs"
        );
    }

    const double ratio =
        inlet_relative_velocity /
        outlet_relative_velocity;

    if (!std::isfinite(ratio)) {
        throw std::runtime_error(
            "Non-finite diffusion ratio"
        );
    }

    const double excess =
        std::max(
            0.0,
            ratio - params.critical_diffusion_ratio
        );

    const double loss =
        0.5 *
        params.loss_coefficient *
        inlet_relative_velocity *
        inlet_relative_velocity *
        excess *
        excess;

    if (!std::isfinite(loss)) {
        throw std::runtime_error(
            "Non-finite diffusion loss"
        );
    }

    return {
        inlet_relative_velocity,
        outlet_relative_velocity,
        ratio,
        excess,
        loss
    };
}

} // namespace turbocompgpu
