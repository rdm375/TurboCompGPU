#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

// Steady, adiabatic, single-inlet/single-outlet impeller.
//
// All quantities are dimensional SI values.
//
// Mass flow is positive in the nominal inlet-to-outlet direction.
// Shaft torque is positive when the shaft does work on the fluid.

struct ConservationGeometry {
    double inlet_radius;   // m
    double outlet_radius;  // m
};

struct ConservationState {
    double inlet_mass_flow;       // kg/s
    double outlet_mass_flow;      // kg/s

    double inlet_whirl;           // m/s
    double outlet_whirl;          // m/s

    double inlet_stagnation_h;    // J/kg
    double outlet_stagnation_h;   // J/kg

    double shaft_torque;          // N m
    double angular_speed;         // rad/s
};

struct ConservationResiduals {
    double mass;             // kg/s
    double angular_momentum; // N m
    double energy;           // W
};

struct ConservationScales {
    double mass;             // kg/s
    double angular_momentum; // N m
    double energy;           // W
};

struct NormalizedConservationResiduals {
    double mass;
    double angular_momentum;
    double energy;
};

inline void validate_conservation_inputs(
    const ConservationGeometry& geometry,
    const ConservationState& state)
{
    const double values[] = {
        geometry.inlet_radius,
        geometry.outlet_radius,

        state.inlet_mass_flow,
        state.outlet_mass_flow,

        state.inlet_whirl,
        state.outlet_whirl,

        state.inlet_stagnation_h,
        state.outlet_stagnation_h,

        state.shaft_torque,
        state.angular_speed
    };

    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "Non-finite conservation input");
        }
    }

    if (geometry.inlet_radius < 0.0 ||
        geometry.outlet_radius <= geometry.inlet_radius) {
        throw std::invalid_argument(
            "Invalid impeller radii");
    }

    if (state.inlet_mass_flow < 0.0 ||
        state.outlet_mass_flow < 0.0) {
        throw std::invalid_argument(
            "Negative mass flow is not supported");
    }

    if (state.angular_speed <= 0.0) {
        throw std::invalid_argument(
            "Angular speed must be positive");
    }
}

inline ConservationResiduals evaluate_conservation(
    const ConservationGeometry& geometry,
    const ConservationState& state)
{
    validate_conservation_inputs(geometry, state);

    const double angular_momentum_in =
        state.inlet_mass_flow *
        geometry.inlet_radius *
        state.inlet_whirl;

    const double angular_momentum_out =
        state.outlet_mass_flow *
        geometry.outlet_radius *
        state.outlet_whirl;

    const double energy_in =
        state.inlet_mass_flow *
        state.inlet_stagnation_h;

    const double energy_out =
        state.outlet_mass_flow *
        state.outlet_stagnation_h;

    const double shaft_power =
        state.shaft_torque * state.angular_speed;

    ConservationResiduals result{
        state.outlet_mass_flow -
            state.inlet_mass_flow,

        state.shaft_torque -
            (angular_momentum_out -
             angular_momentum_in),

        energy_out -
            energy_in -
            shaft_power
    };

    if (!std::isfinite(result.mass) ||
        !std::isfinite(result.angular_momentum) ||
        !std::isfinite(result.energy)) {
        throw std::overflow_error(
            "Non-finite conservation residual");
    }

    return result;
}

// Convert independently specified physical residual scales
// into dimensionless residuals.
//
// These scales are intended to be fixed for a nonlinear solve,
// not recomputed from every trial state.

inline NormalizedConservationResiduals
normalize_conservation(
    const ConservationResiduals& residuals,
    const ConservationScales& scales)
{
    if (!std::isfinite(scales.mass) ||
        !std::isfinite(scales.angular_momentum) ||
        !std::isfinite(scales.energy) ||
        scales.mass <= 0.0 ||
        scales.angular_momentum <= 0.0 ||
        scales.energy <= 0.0) {
        throw std::invalid_argument(
            "Invalid conservation scales");
    }

    const NormalizedConservationResiduals result{
        residuals.mass / scales.mass,
        residuals.angular_momentum /
            scales.angular_momentum,
        residuals.energy / scales.energy
    };

    if (!std::isfinite(result.mass) ||
        !std::isfinite(result.angular_momentum) ||
        !std::isfinite(result.energy)) {
        throw std::overflow_error(
            "Non-finite normalized conservation residual");
    }

    return result;
}

// Construct fixed residual scales from nominal operating values.
//
// The specific-work scale is an independent user-supplied
// characteristic value, not a residual-derived quantity.

inline ConservationScales make_conservation_scales(
    double nominal_mass_flow,
    double nominal_specific_work,
    double nominal_angular_speed)
{
    if (!std::isfinite(nominal_mass_flow) ||
        !std::isfinite(nominal_specific_work) ||
        !std::isfinite(nominal_angular_speed) ||
        nominal_mass_flow <= 0.0 ||
        nominal_specific_work <= 0.0 ||
        nominal_angular_speed <= 0.0) {
        throw std::invalid_argument(
            "Invalid nominal conservation scales");
    }

    const double nominal_power =
        nominal_mass_flow * nominal_specific_work;

    const double nominal_torque =
        nominal_power / nominal_angular_speed;

    if (!std::isfinite(nominal_power) ||
        !std::isfinite(nominal_torque) ||
        nominal_torque <= 0.0) {
        throw std::overflow_error(
            "Invalid derived conservation scales");
    }

    return {
        nominal_mass_flow,
        nominal_torque,
        nominal_power
    };
}

} // namespace turbocompgpu
