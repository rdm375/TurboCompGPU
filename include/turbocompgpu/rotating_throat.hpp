#pragma once

#include <turbocompgpu/choking_mass_flux.hpp>

#include <thermomodel/model.hpp>
#include <thermomodel/derivatives.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

// Locally one-dimensional, inviscid, adiabatic rotating throat.
//
// W_n: velocity normal to the throat.
// W_t: prescribed tangential relative velocity.
// U:   blade speed at the throat.
//
// Rothalpy:
//     I = h + (W_n^2 + W_t^2)/2 - U^2/2
//
// Critical condition:
//     W_n = a
//
// Consequently:
//     h + a^2/2 = I + (U^2 - W_t^2)/2
//
// The supplied effective stagnation state defines:
//
//     h0_eff = I + (U^2 - W_t^2)/2
//     s0_eff = s
//
// This representation permits direct reuse of the validated
// stationary choking solver.
//
// Assumptions:
//   - steady, adiabatic, isentropic flow
//   - constant U and W_t across the local throat calculation
//   - normal sonic choking
//   - prescribed effective throat area
//
// No compressor-wide choke boundary is inferred here.

struct RotatingThroatInputs {
    double effective_stagnation_temperature;
    double effective_stagnation_pressure;

    double radius;
    double rotational_speed_rpm;
    double tangential_relative_velocity;
    double effective_area;
};

struct RotatingThroatResult {
    ChokingMassFluxResult critical;

    double blade_speed;
    double rothalpy;
    double normal_relative_velocity;
    double tangential_relative_velocity;
    double total_relative_velocity;

    double critical_mass_flow;
    double rotating_energy_residual;
};

inline RotatingThroatResult evaluate_rotating_throat(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const RotatingThroatInputs& input)
{
    const auto valid = [](double x) {
        return std::isfinite(x);
    };

    if (!(valid(input.radius) &&
          input.radius >= 0.0 &&
          valid(input.rotational_speed_rpm) &&
          input.rotational_speed_rpm >= 0.0 &&
          valid(input.tangential_relative_velocity) &&
          valid(input.effective_area) &&
          input.effective_area > 0.0)) {
        throw std::invalid_argument(
            "Invalid rotating throat geometry or velocity");
    }

    const double omega =
        input.rotational_speed_rpm *
        (2.0 * std::numbers::pi / 60.0);

    const double U = omega * input.radius;
    const double Wt = input.tangential_relative_velocity;

    if (!valid(U)) {
        throw std::invalid_argument(
            "Invalid rotating throat blade speed");
    }

    const auto stagnation = model.evaluate_Tp(
        input.effective_stagnation_temperature,
        input.effective_stagnation_pressure);

    if (!stagnation) {
        throw std::runtime_error(
            "Effective stagnation state evaluation failed");
    }

    const double h0 = stagnation.value.h;

    const double I =
        h0 - 0.5 * (U * U - Wt * Wt);

    if (!valid(I)) {
        throw std::runtime_error(
            "Nonfinite rotating throat rothalpy");
    }

    const auto critical = evaluate_choking_mass_flux(
        model,
        derivatives,
        input.effective_stagnation_temperature,
        input.effective_stagnation_pressure);

    const double Wn = critical.sound_speed;

    const double W2 = Wn * Wn + Wt * Wt;

    const double rotating_energy_residual =
        critical.enthalpy +
        0.5 * W2 -
        0.5 * U * U -
        I;

    return {
        critical,
        U,
        I,
        Wn,
        Wt,
        std::sqrt(W2),
        input.effective_area * critical.critical_mass_flux,
        rotating_energy_residual
    };
}

} // namespace turbocompgpu
