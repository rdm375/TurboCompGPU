#pragma once

#include <thermomodel/model.hpp>

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct EntropyLossModel {
    double entropy_generation; // J/(kg K)
};

struct EntropyLossResult {
    double inlet_stagnation_entropy;
    double outlet_stagnation_entropy;

    double actual_entropy_increase;
    double prescribed_entropy_increase;

    double entropy_residual; // J/(kg K)
};

inline EntropyLossResult evaluate_entropy_loss(
    const thermomodel::State& inlet_stagnation,
    const thermomodel::State& outlet_stagnation,
    const EntropyLossModel& loss)
{
    if (!std::isfinite(loss.entropy_generation) ||
        loss.entropy_generation < 0.0) {
        throw std::invalid_argument(
            "Entropy generation must be finite and nonnegative");
    }

    if (!std::isfinite(inlet_stagnation.s) ||
        !std::isfinite(outlet_stagnation.s)) {
        throw std::invalid_argument(
            "Non-finite stagnation entropy");
    }

    const double actual =
        outlet_stagnation.s -
        inlet_stagnation.s;

    const double residual =
        actual - loss.entropy_generation;

    if (!std::isfinite(actual) ||
        !std::isfinite(residual)) {
        throw std::runtime_error(
            "Non-finite entropy-loss result");
    }

    return {
        inlet_stagnation.s,
        outlet_stagnation.s,
        actual,
        loss.entropy_generation,
        residual
    };
}

} // namespace turbocompgpu
