#pragma once

#include <turbocompgpu/thermodynamic_solver.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

struct AbsoluteVelocity {
    double meridional = 0.0;  // m/s
    double tangential = 0.0;  // m/s
};

struct FlowStation {
    thermomodel::State static_state;
    AbsoluteVelocity velocity;
};

struct StagnationResult {
    thermomodel::State stagnation_state;

    double stagnation_enthalpy = 0.0; // J/kg
    double stagnation_entropy = 0.0;  // J/(kg K)

    double enthalpy_residual = 0.0;
    double entropy_residual = 0.0;

    int pressure_iterations = 0;
    int pressure_evaluations = 0;
    int temperature_evaluations = 0;

    bool resolution_limited = false;
};

inline double absolute_speed_squared(
    const AbsoluteVelocity& velocity)
{
    if (!std::isfinite(velocity.meridional) ||
        !std::isfinite(velocity.tangential)) {
        throw std::invalid_argument(
            "Non-finite absolute velocity");
    }

    const double V = std::hypot(
        velocity.meridional,
        velocity.tangential);

    const double V2 = V * V;

    if (!std::isfinite(V2)) {
        throw std::overflow_error(
            "Absolute velocity squared overflow");
    }

    return V2;
}

inline double stagnation_enthalpy(
    const FlowStation& station)
{
    if (!std::isfinite(station.static_state.h)) {
        throw std::invalid_argument(
            "Non-finite static enthalpy");
    }

    const double h0 =
        station.static_state.h +
        0.5 * absolute_speed_squared(station.velocity);

    if (!std::isfinite(h0)) {
        throw std::overflow_error(
            "Stagnation enthalpy overflow");
    }

    return h0;
}

inline StagnationResult evaluate_stagnation_state(
    const thermomodel::Model& model,
    const FlowStation& station)
{
    const auto& supplied = station.static_state;

    // Re-evaluate the static state using the supplied model.
    // This prevents inconsistent cached thermodynamic properties
    // from silently entering the stagnation calculation.
    const auto static_state =
        evaluate(model, supplied.T, supplied.p);

    const double V2 =
        absolute_speed_squared(station.velocity);

    const double h0 =
        static_state.h + 0.5 * V2;

    if (!std::isfinite(h0)) {
        throw std::overflow_error(
            "Stagnation enthalpy overflow");
    }

    const double s0 = static_state.s;

    StagnationResult result{};
    result.stagnation_enthalpy = h0;
    result.stagnation_entropy = s0;

    // Zero velocity: stagnation and static states coincide.
    if (V2 == 0.0) {
        result.stagnation_state = static_state;
        return result;
    }

    const auto limits = model.metadata();

    if (!(limits.pmax >= static_state.p) ||
        !std::isfinite(limits.pmax)) {
        throw std::runtime_error(
            "Invalid stagnation pressure search limits");
    }

    struct Trial {
        thermomodel::State state;
        double entropy_residual;
        bool resolution_limited;
    };

    auto trial_at_pressure = [&](double pressure) -> Trial {
        // The fixed-pressure enthalpy inversion is our
        // previously validated safeguarded Newton solver.
        const auto inversion = invert_temperature(
            model,
            pressure,
            h0,
            Property::Enthalpy,
            static_state.T);

        ++result.pressure_evaluations;
        result.temperature_evaluations +=
            inversion.evaluations;

        return {
            inversion.state,
            inversion.state.s - s0,
            inversion.resolution_limited
        };
    };

    // Use logarithmic pressure relative to the static state.
    const double x_max =
        std::log(limits.pmax / static_state.p);

    if (!(x_max > 0.0)) {
        throw std::runtime_error(
            "No pressure range available for stagnation inversion");
    }

    // At static pressure, adding kinetic energy increases
    // temperature and entropy on the stable vapor branch.
    const Trial initial =
        trial_at_pressure(static_state.p);

    if (!(initial.entropy_residual > 0.0)) {
        throw std::runtime_error(
            "Unexpected entropy ordering at static pressure");
    }

    double x_lo = 0.0;
    double f_lo = initial.entropy_residual;

    double x_hi = 0.0;
    double f_hi = f_lo;

    Trial upper = initial;
    bool bracketed = false;

    double step = std::min(0.02, x_max);

    // Geometrically expand the logarithmic pressure interval.
    for (int expansion = 0; expansion < 80; ++expansion) {
        const double candidate =
            std::min(x_max, x_hi + step);

        if (candidate <= x_hi) {
            break;
        }

        const double pressure =
            static_state.p * std::exp(candidate);

        // An invalid intermediate thermodynamic state is
        // reported, not silently skipped across a phase boundary.
        upper = trial_at_pressure(pressure);

        x_hi = candidate;
        f_hi = upper.entropy_residual;

        if (f_hi <= 0.0) {
            bracketed = true;
            break;
        }

        step *= 2.0;
    }

    if (!bracketed) {
        throw std::runtime_error(
            "Stagnation pressure root not bracketed");
    }

    // The entropy residual along an isenthalpic curve obeys
    //
    //     ds/d(ln p)|h = -v*p/T
    //
    // for a simple compressible substance.
    //
    // Thus an entropy residual corresponds approximately to
    // a logarithmic-pressure correction:
    //
    //     dx = residual * rho*T/p.
    //
    // This is also reference-state invariant.

    constexpr double relative_pressure_tolerance = 1.0e-12;
    constexpr int max_iterations = 200;

    // Begin with a bracketed iterate. Subsequent iterations
    // use safeguarded Newton corrections whenever possible.
    double x = x_lo + 0.5 * (x_hi - x_lo);

    const double initial_newton =
        initial.entropy_residual *
        initial.state.rho *
        initial.state.T /
        static_state.p;

    if (std::isfinite(initial_newton) &&
        initial_newton > x_lo + 0.05 * (x_hi - x_lo) &&
        initial_newton < x_hi - 0.05 * (x_hi - x_lo)) {
        x = initial_newton;
    }

    for (int iteration = 1;
         iteration <= max_iterations;
         ++iteration) {

        result.pressure_iterations = iteration;

        const double pressure =
            static_state.p * std::exp(x);

        const Trial current =
            trial_at_pressure(pressure);

        const double f = current.entropy_residual;

        if (!std::isfinite(f)) {
            throw std::runtime_error(
                "Non-finite stagnation entropy residual");
        }

        // Exact derivative along the constant-enthalpy curve:
        //
        // dF/dx = -p/(rho*T), x = ln(p/p_static).
        //
        // The corresponding Newton correction is
        // dx = F*rho*T/p.
        const double correction =
            f * current.state.rho *
            current.state.T / pressure;

        if (!std::isfinite(correction)) {
            throw std::runtime_error(
                "Non-finite stagnation pressure correction");
        }

        const double pressure_equivalent_residual =
            std::abs(correction);

        if (pressure_equivalent_residual <=
            relative_pressure_tolerance) {

            result.stagnation_state = current.state;

            result.enthalpy_residual =
                current.state.h - h0;

            result.entropy_residual = f;

            result.resolution_limited =
                current.resolution_limited;

            return result;
        }

        // Preserve the sign-changing bracket.
        if (f > 0.0) {
            x_lo = x;
            f_lo = f;
        } else {
            x_hi = x;
            f_hi = f;
        }

        const double width = x_hi - x_lo;

        // Candidate Newton iterate.
        const double x_newton = x + correction;

        // Accept any representable Newton step strictly
        // inside the bracket. Otherwise use bisection.
        if (std::isfinite(x_newton) &&
            x_newton > x_lo &&
            x_newton < x_hi &&
            x_newton != x) {

            x = x_newton;

        } else {
            x = x_lo + 0.5 * width;
        }

        // Do not accept a small bracket as convergence.
        // Detect inability to represent a new iterate.
        if (x <= x_lo || x >= x_hi) {
            throw std::runtime_error(
                "Stagnation pressure inversion stagnated "
                "before meeting the residual criterion");
        }
    }

    throw std::runtime_error(
        "Stagnation pressure inversion did not converge");
}

} // namespace turbocompgpu
