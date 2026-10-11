#pragma once

#include <thermomodel/model.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

enum class Property {
    Enthalpy,
    Entropy
};

struct InversionResult {
    thermomodel::State state;
    int iterations = 0;
    int evaluations = 0;
    double residual = 0.0;
    double bracket_width = 0.0;
    bool resolution_limited = false;
};

inline bool vapor_like(thermomodel::Phase phase)
{
    return phase == thermomodel::Phase::Vapor ||
           phase == thermomodel::Phase::Supercritical;
}

inline thermomodel::State evaluate(
    const thermomodel::Model& model,
    double T,
    double pressure)
{
    if (!std::isfinite(T) ||
        !std::isfinite(pressure) ||
        T <= 0.0 ||
        pressure <= 0.0) {
        throw std::invalid_argument("Invalid thermodynamic inputs");
    }

    const auto limits = model.metadata();

    if (T < limits.Tmin ||
        T > limits.Tmax ||
        pressure < limits.pmin ||
        pressure > limits.pmax) {
        throw std::out_of_range(
            "Thermodynamic state outside model limits");
    }

    const auto result = model.evaluate_Tp(T, pressure);

    if (!result) {
        throw std::runtime_error(
            "Thermodynamic evaluation failed: " + result.message);
    }

    const auto& state = result.value;

    if (!std::isfinite(state.T) ||
        !std::isfinite(state.p) ||
        !std::isfinite(state.h) ||
        !std::isfinite(state.s) ||
        !std::isfinite(state.rho) ||
        state.rho <= 0.0) {
        throw std::runtime_error(
            "Invalid thermodynamic state");
    }

    if (!vapor_like(state.phase)) {
        throw std::runtime_error(
            "Thermodynamic state is not vapor-like");
    }

    return state;
}

inline double property_value(
    const thermomodel::State& state,
    Property property)
{
    return property == Property::Enthalpy
        ? state.h
        : state.s;
}

inline InversionResult invert_temperature(
    const thermomodel::Model& model,
    double pressure,
    double target,
    Property property,
    double initial_temperature)
{
    if (!std::isfinite(target) ||
        !std::isfinite(initial_temperature)) {
        throw std::invalid_argument(
            "Invalid temperature inversion inputs");
    }

    const auto limits = model.metadata();

    if (!(limits.Tmin > 0.0) ||
        !(limits.Tmax > limits.Tmin)) {
        throw std::runtime_error(
            "Invalid thermodynamic model temperature limits");
    }

    int evaluations = 0;

    auto evaluate_counted = [&](double T) {
        auto state = evaluate(model, T, pressure);
        ++evaluations;
        return state;
    };

    const auto initial =
        evaluate_counted(initial_temperature);

    const double initial_value =
        property_value(initial, property);

    const double initial_residual =
        initial_value - target;

    constexpr double relative_temperature_tolerance = 1.0e-13;
    constexpr int max_iterations = 200;
    constexpr int max_expansions = 80;

    // An exact initial match requires no iteration.
    if (initial_residual == 0.0) {
        return {
            initial, 0, evaluations,
            initial_residual, 0.0
        };
    }

    // The current implementation assumes that enthalpy and
    // entropy increase monotonically with temperature
    // along the selected single-phase branch.
    const bool increase_temperature =
        initial_residual < 0.0;

    double lo = initial_temperature;
    double hi = initial_temperature;

    double f_lo = initial_residual;
    double f_hi = initial_residual;

    double step = std::max(
        1.0, 0.02 * initial_temperature);

    bool bracketed = false;

    for (int expansion = 0;
         expansion < max_expansions;
         ++expansion) {

        const double candidate = increase_temperature
            ? std::min(limits.Tmax, hi + step)
            : std::max(limits.Tmin, lo - step);

        const double previous = increase_temperature
            ? hi : lo;

        if (candidate == previous) {
            break;
        }

        thermomodel::State state;

        try {
            state = evaluate_counted(candidate);
        }
        catch (const std::exception& error) {
            throw std::runtime_error(
                std::string(
                    "Temperature bracketing encountered "
                    "an invalid state: ") + error.what());
        }

        const double residual =
            property_value(state, property) - target;

        if (increase_temperature) {
            hi = candidate;
            f_hi = residual;
        } else {
            lo = candidate;
            f_lo = residual;
        }

        if (residual == 0.0) {
            return {
                state, 0, evaluations,
                residual, hi - lo
            };
        }

        if ((f_lo < 0.0) != (f_hi < 0.0)) {
            bracketed = true;
            break;
        }

        step *= 2.0;
    }

    if (!bracketed) {
        throw std::runtime_error(
            "Temperature inversion target not bracketed");
    }


    // Safeguarded Newton iteration.
    //
    // The bracket is retained throughout. Newton steps are
    // accepted only when they remain safely inside it.

    auto derivative = [&](const thermomodel::State& state) {
        if (!std::isfinite(state.cp) || state.cp <= 0.0) {
            throw std::runtime_error(
                "Invalid heat capacity during temperature inversion");
        }

        return property == Property::Enthalpy
            ? state.cp
            : state.cp / state.T;
    };

    double T = lo + 0.5 * (hi - lo);

    for (int iteration = 1;
         iteration <= max_iterations;
         ++iteration) {

        const auto state = evaluate_counted(T);

        const double residual =
            property_value(state, property) - target;

        const double slope = derivative(state);

        if (!std::isfinite(slope) || slope <= 0.0) {
            throw std::runtime_error(
                "Invalid thermodynamic property derivative");
        }

        const double temperature_residual =
            std::abs(residual) / slope;

        const double temperature_tolerance =
            relative_temperature_tolerance *
            std::max(1.0, std::abs(T));

        // Temperature-equivalent convergence is invariant under
        // exact additive shifts of the property reference.
        //
        // The roundoff allowance recognizes that subtraction of
        // large absolute property values can limit resolution.
        constexpr double roundoff_factor = 8.0;

        const double roundoff_allowance =
            roundoff_factor *
            std::numeric_limits<double>::epsilon() *
            std::max({
                1.0,
                std::abs(property_value(state, property)),
                std::abs(target)
            });

        const bool resolution_limited =
            temperature_residual > temperature_tolerance &&
            std::abs(residual) <= roundoff_allowance;

        if (temperature_residual <= temperature_tolerance ||
            resolution_limited) {
            return {
                state,
                iteration,
                evaluations,
                residual,
                hi - lo,
                resolution_limited
            };
        }

        if ((residual < 0.0) == (f_lo < 0.0)) {
            lo = T;
            f_lo = residual;
        } else {
            hi = T;
            f_hi = residual;
        }

        const double width = hi - lo;

        if (width <= temperature_tolerance) {
            throw std::runtime_error(
                "Temperature inversion reached floating-point "
                "resolution before residual convergence");
        }

        const double newton_step = -residual / slope;
        const double newton_candidate = T + newton_step;

        // Accept a finite Newton step strictly inside the bracket.
        // Otherwise fall back to bisection.
        if (std::isfinite(newton_candidate) &&
            newton_candidate > lo &&
            newton_candidate < hi &&
            newton_candidate != T) {

            T = newton_candidate;

        } else {

            T = lo + 0.5 * width;
        }

        if (T == lo || T == hi) {
            throw std::runtime_error(
                "Temperature inversion stagnated");
        }
    }

    throw std::runtime_error(
        "Temperature inversion exceeded iteration limit");
}

} // namespace turbocompgpu
