#pragma once

#include <turbocompgpu/station.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

// Stage 8B.1: adiabatic, lossless stationary-component reference.
//
// The outlet static pressure is prescribed. The resulting outlet
// speed and required effective flow area are calculated, not imposed.
//
// This is a thermodynamic reference component, not a predictive
// diffuser or collector geometry model.

enum class StationaryComponentKind {
    Diffuser,
    Collector
};

struct StageStation {
    thermomodel::State static_state;
    StagnationResult stagnation;
    double absolute_speed;
    double mass_flow;
};

struct StationaryComponentResult {
    StageStation inlet;
    StageStation outlet;

    double required_effective_area;

    double mass_residual;
    double stagnation_enthalpy_residual;
    double entropy_residual;
    double stagnation_pressure_residual;
};

inline StageStation make_stage_station(
    const thermomodel::Model& model,
    const thermomodel::State& static_state,
    double absolute_speed,
    double mass_flow)
{
    if (!std::isfinite(absolute_speed) ||
        absolute_speed < 0.0 ||
        !std::isfinite(mass_flow) ||
        mass_flow <= 0.0) {
        throw std::invalid_argument(
            "Invalid stage station speed or mass flow");
    }

    const auto state = evaluate(
        model, static_state.T, static_state.p);

    const FlowStation flow{
        state, {absolute_speed, 0.0}
    };

    return {
        state,
        evaluate_stagnation_state(model, flow),
        absolute_speed,
        mass_flow
    };
}

inline StationaryComponentResult
evaluate_lossless_stationary_component(
    const thermomodel::Model& model,
    const StageStation& inlet,
    double outlet_static_pressure,
    StationaryComponentKind kind)
{
    if (!std::isfinite(outlet_static_pressure) ||
        outlet_static_pressure <= 0.0) {
        throw std::invalid_argument(
            "Invalid stationary-component outlet pressure");
    }

    // This reference is restricted to pressure recovery.
    if (outlet_static_pressure < inlet.static_state.p) {
        throw std::invalid_argument(
            "Stationary reference requires nondecreasing static pressure");
    }

    const double h0 =
        inlet.stagnation.stagnation_enthalpy;

    const double s0 =
        inlet.static_state.s;

    const double p0 =
        inlet.stagnation.stagnation_state.p;

    if (outlet_static_pressure > p0) {
        throw std::invalid_argument(
            "Outlet static pressure exceeds inlet stagnation pressure");
    }

    // Isentropic static state at the prescribed pressure.
    const auto inversion = invert_temperature(
        model,
        outlet_static_pressure,
        s0,
        Property::Entropy,
        inlet.static_state.T
    );

    const auto& outlet_static = inversion.state;

    const double kinetic_energy =
        h0 - outlet_static.h;

    // Allow only roundoff-scale negative kinetic energy.
    const double kinetic_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max({
            1.0,
            std::abs(h0),
            std::abs(outlet_static.h)
        });

    if (kinetic_energy < -kinetic_tolerance) {
        throw std::runtime_error(
            "Negative outlet kinetic energy");
    }

    const double outlet_speed =
        std::sqrt(2.0 * std::max(0.0, kinetic_energy));

    // This is a diffuser/collector reference: no acceleration.
    const double speed_tolerance =
        1.0e-9 * std::max(1.0, inlet.absolute_speed);

    if (outlet_speed > inlet.absolute_speed + speed_tolerance) {
        throw std::runtime_error(
            "Stationary reference unexpectedly accelerates flow");
    }

    if (!(outlet_speed > 0.0)) {
        throw std::runtime_error(
            "Zero outlet speed requires infinite flow area");
    }

    const double area =
        inlet.mass_flow /
        (outlet_static.rho * outlet_speed);

    if (!std::isfinite(area) || area <= 0.0) {
        throw std::runtime_error(
            "Invalid required effective flow area");
    }

    const auto outlet = make_stage_station(
        model,
        outlet_static,
        outlet_speed,
        inlet.mass_flow
    );

    return {
        inlet,
        outlet,
        area,
        outlet.mass_flow - inlet.mass_flow,
        outlet.stagnation.stagnation_enthalpy -
            inlet.stagnation.stagnation_enthalpy,
        outlet.static_state.s - inlet.static_state.s,
        outlet.stagnation.stagnation_state.p -
            inlet.stagnation.stagnation_state.p
    };
}

struct CompleteStageReference {
    StageStation inlet;
    StageStation impeller_exit;

    StationaryComponentResult diffuser;
    StationaryComponentResult collector;

    double fluid_specific_work;
    double fluid_power;
    double total_to_total_pressure_ratio;

    double mass_residual;
    double energy_residual;
    double stationary_entropy_change;
};

inline CompleteStageReference evaluate_complete_stage_reference(
    const thermomodel::Model& model,
    const StageStation& inlet,
    const StageStation& impeller_exit,
    double diffuser_exit_static_pressure,
    double collector_exit_static_pressure)
{
    if (!(inlet.mass_flow > 0.0) ||
        !(impeller_exit.mass_flow > 0.0)) {
        throw std::invalid_argument(
            "Stage mass flow must be positive");
    }

    const double mass_residual =
        impeller_exit.mass_flow - inlet.mass_flow;

    const double mass_scale =
        std::max(inlet.mass_flow, impeller_exit.mass_flow);

    if (std::abs(mass_residual) >
        1.0e-12 * mass_scale) {
        throw std::invalid_argument(
            "Impeller inlet and exit mass flows differ");
    }

    const auto diffuser =
        evaluate_lossless_stationary_component(
            model,
            impeller_exit,
            diffuser_exit_static_pressure,
            StationaryComponentKind::Diffuser
        );

    const auto collector =
        evaluate_lossless_stationary_component(
            model,
            diffuser.outlet,
            collector_exit_static_pressure,
            StationaryComponentKind::Collector
        );

    const double work =
        impeller_exit.stagnation.stagnation_enthalpy -
        inlet.stagnation.stagnation_enthalpy;

    const double final_work =
        collector.outlet.stagnation.stagnation_enthalpy -
        inlet.stagnation.stagnation_enthalpy;

    return {
        inlet,
        impeller_exit,
        diffuser,
        collector,
        work,
        inlet.mass_flow * work,
        collector.outlet.stagnation.stagnation_state.p /
            inlet.stagnation.stagnation_state.p,
        collector.outlet.mass_flow - inlet.mass_flow,
        final_work - work,
        collector.outlet.static_state.s -
            impeller_exit.static_state.s
    };
}

} // namespace turbocompgpu
