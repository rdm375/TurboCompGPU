#pragma once

#include <turbocompgpu/isentropic_sound_speed.hpp>

#include <thermomodel/model.hpp>
#include <thermomodel/derivatives.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Stationary, inviscid, adiabatic, isentropic throat.
//
// h0 = h + V^2/2
// s  = s0
// V  = a at the critical state
// G* = rho* a*
//
// This is a standalone thermodynamic kernel, not yet
// a rotating-impeller choking model.

struct ChokingMassFluxResult {
    double temperature;
    double pressure;
    double density;
    double enthalpy;
    double entropy;
    double sound_speed;
    double critical_mass_flux;
    double energy_residual;
    double entropy_residual;
    int iterations;
};

inline ChokingMassFluxResult evaluate_choking_mass_flux(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    double stagnation_temperature,
    double stagnation_pressure)
{
    if (!(std::isfinite(stagnation_temperature) &&
          stagnation_temperature > 0.0 &&
          std::isfinite(stagnation_pressure) &&
          stagnation_pressure > 0.0)) {
        throw std::invalid_argument(
            "Invalid choking calculation inputs");
    }

    const auto stagnation =
        model.evaluate_Tp(
            stagnation_temperature,
            stagnation_pressure);

    if (!stagnation) {
        throw std::runtime_error(
            "Choking stagnation state evaluation failed");
    }

    const double h0 = stagnation.value.h;
    const double s0 = stagnation.value.s;

    const auto limits = model.metadata();

    // For a trial pressure, solve s(T,p)=s0 using
    // safeguarded bisection in temperature.
    const auto isentropic_state =
        [&](double pressure) -> thermomodel::State {

        double lo = limits.Tmin;
        double hi = limits.Tmax;

        auto slo = model.evaluate_Tp(lo, pressure);
        auto shi = model.evaluate_Tp(hi, pressure);

        if (!slo || !shi) {
            throw std::runtime_error(
                "Isentropic temperature bracket invalid");
        }

        double flo = slo.value.s - s0;
        double fhi = shi.value.s - s0;

        if (!(flo <= 0.0 && fhi >= 0.0)) {
            throw std::runtime_error(
                "Isentropic temperature root not bracketed");
        }

        for (int iteration = 0; iteration < 120; ++iteration) {
            const double mid = 0.5 * (lo + hi);

            const auto state =
                model.evaluate_Tp(mid, pressure);

            if (!state) {
                throw std::runtime_error(
                    "Isentropic state evaluation failed");
            }

            const double f = state.value.s - s0;

            if (std::abs(f) <= 1.0e-10 ||
                (hi - lo) <=
                    1.0e-12 * std::max(1.0, mid)) {
                return state.value;
            }

            if (f < 0.0) {
                lo = mid;
                flo = f;
            } else {
                hi = mid;
                fhi = f;
            }
        }

        throw std::runtime_error(
            "Isentropic temperature solve did not converge");
    };

    // Sonic energy residual:
    //
    // F(p) = h(p,s0) + a(p,s0)^2/2 - h0.
    //
    // At p0, F>0. At sufficiently low p, F<0.

    const auto residual =
        [&](double pressure) {

        const auto state = isentropic_state(pressure);

        const auto sound =
            evaluate_isentropic_sound_speed(
                derivatives, state.T, pressure);

        return state.h +
               0.5 * sound.sound_speed_squared - h0;
    };

    double upper = stagnation_pressure;
    double lower = 0.5 * upper;

    const double minimum_pressure =
        std::max(limits.pmin, 1.0);

    while (lower > minimum_pressure &&
           residual(lower) > 0.0) {
        lower *= 0.5;
    }

    lower = std::max(lower, minimum_pressure);

    double flo = residual(lower);
    double fhi = residual(upper);

    if (!(flo <= 0.0 && fhi >= 0.0)) {
        throw std::runtime_error(
            "Sonic pressure root not bracketed");
    }

    for (int iteration = 1; iteration <= 120; ++iteration) {
        const double pressure = 0.5 * (lower + upper);

        const auto state = isentropic_state(pressure);

        const auto sound =
            evaluate_isentropic_sound_speed(
                derivatives, state.T, pressure);

        const double f =
            state.h +
            0.5 * sound.sound_speed_squared - h0;

        const double energy_scale =
            std::max({
                1.0,
                std::abs(h0),
                std::abs(state.h),
                0.5 * sound.sound_speed_squared
            });

        // Terminate on the requested energy accuracy,
        // or when the pressure bracket is exhausted.
        //
        // The returned energy_residual remains the
        // independently computed physical residual.
        const double pressure_scale =
            std::max(1.0, std::abs(pressure));

        const bool energy_converged =
            std::abs(f) <= 1.0e-12 * energy_scale;

        const bool pressure_converged =
            (upper - lower) <=
            1.0e-13 * pressure_scale;

        if (energy_converged || pressure_converged) {

            return {
                state.T,
                pressure,
                state.rho,
                state.h,
                state.s,
                sound.sound_speed,
                state.rho * sound.sound_speed,
                f,
                state.s - s0,
                iteration
            };
        }

        if (f < 0.0) {
            lower = pressure;
            flo = f;
        } else {
            upper = pressure;
            fhi = f;
        }
    }

    throw std::runtime_error(
        "Sonic pressure solve did not converge");
}

} // namespace turbocompgpu
