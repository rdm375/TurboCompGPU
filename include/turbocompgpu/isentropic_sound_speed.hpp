#pragma once

#include <thermomodel/derivatives.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

struct IsentropicSoundSpeedResult {
    double sound_speed;
    double sound_speed_squared;
    double isentropic_density_derivative;
};

// Calculate a^2 = (dp/drho)_s using analytical EOS derivatives.
//
// Along an isentrope:
//
//   dT/dp|s = -(ds/dp|T) / (ds/dT|p)
//
// Therefore:
//
//   drho/dp|s = drho/dp|T
//              - drho/dT|p * (ds/dp|T)/(ds/dT|p)
//
//   a^2 = 1 / (drho/dp|s)
//
// The calculation is independent of enthalpy and entropy
// reference-state offsets.

inline IsentropicSoundSpeedResult
evaluate_isentropic_sound_speed(
    const thermomodel::DifferentiableModel& model,
    double temperature,
    double pressure)
{
    if (!(std::isfinite(temperature) && temperature > 0.0 &&
          std::isfinite(pressure) && pressure > 0.0)) {
        throw std::invalid_argument(
            "Invalid thermodynamic state for sound speed");
    }

    const auto result =
        model.derivatives_Tp(temperature, pressure);

    if (!result) {
        throw std::runtime_error(
            "Thermodynamic derivative evaluation failed: " +
            result.message);
    }

    const auto& d = result.value;

    if (!(std::isfinite(d.drho_dp_T) &&
          std::isfinite(d.drho_dT_p) &&
          std::isfinite(d.ds_dp_T) &&
          std::isfinite(d.ds_dT_p) &&
          d.ds_dT_p > 0.0)) {
        throw std::runtime_error(
            "Invalid thermodynamic derivatives for sound speed");
    }

    const double drho_dp_s =
        d.drho_dp_T -
        d.drho_dT_p * (d.ds_dp_T / d.ds_dT_p);

    if (!(std::isfinite(drho_dp_s) && drho_dp_s > 0.0)) {
        throw std::runtime_error(
            "Nonpositive isentropic compressibility");
    }

    const double a2 = 1.0 / drho_dp_s;
    const double a = std::sqrt(a2);

    if (!(std::isfinite(a2) && a2 > 0.0 &&
          std::isfinite(a) && a > 0.0)) {
        throw std::runtime_error(
            "Invalid calculated sound speed");
    }

    return {a, a2, drho_dp_s};
}

} // namespace turbocompgpu
