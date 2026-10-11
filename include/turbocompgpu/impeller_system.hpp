#pragma once

#include <turbocompgpu/continuity.hpp>
#include <turbocompgpu/conservation.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/model.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct ImpellerSystemGeometry {
    InletAnnulus inlet;
    ImpellerExitPassage exit;
    ImpellerGeometry impeller;
};

struct ImpellerTrialState {
    double inlet_temperature;
    double inlet_pressure;

    double outlet_temperature;
    double outlet_pressure;

    double inlet_mass_flow;
    double outlet_mass_flow;

    double rotational_speed_rpm;
    double slip_factor;

    double inlet_prewhirl;
    double outlet_whirl;

    double shaft_torque;
};

struct ImpellerSystemResult {
    thermomodel::State inlet_static;
    thermomodel::State outlet_static;

    ContinuityResult continuity;
    ImpellerResult impeller;

    StagnationResult inlet_stagnation;
    StagnationResult outlet_stagnation;

    ConservationResiduals conservation;

    double predicted_outlet_whirl;
    double trial_outlet_whirl;
    double slip_residual;

    double angular_speed;
    double euler_specific_work;
    double stagnation_enthalpy_change;
};

inline ImpellerSystemResult evaluate_impeller_system(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const ImpellerTrialState& trial)
{
    // Check the geometric relationship between the
    // continuity passage and impeller exit.
    if (!std::isfinite(geometry.exit.diameter) ||
        !std::isfinite(geometry.impeller.exit_diameter) ||
        geometry.exit.diameter !=
            geometry.impeller.exit_diameter) {
        throw std::invalid_argument(
            "Impeller and exit passage diameters differ");
    }

    // The impeller inlet diameter is the representative
    // radius used for inlet angular momentum. It is not
    // necessarily the inlet annulus shroud diameter.

    if (!std::isfinite(trial.outlet_whirl) ||
        !std::isfinite(trial.shaft_torque) ||
        !std::isfinite(trial.inlet_prewhirl)) {
        throw std::invalid_argument(
            "Non-finite independent trial variable");
    }

    // 1. Independent static thermodynamic states.
    const auto inlet_evaluation =
        model.evaluate_Tp(
            trial.inlet_temperature,
            trial.inlet_pressure);

    if (!inlet_evaluation) {
        throw std::runtime_error(
            "Inlet thermodynamic evaluation failed: " +
            inlet_evaluation.message);
    }

    const auto outlet_evaluation =
        model.evaluate_Tp(
            trial.outlet_temperature,
            trial.outlet_pressure);

    if (!outlet_evaluation) {
        throw std::runtime_error(
            "Outlet thermodynamic evaluation failed: " +
            outlet_evaluation.message);
    }

    const auto inlet_static = inlet_evaluation.value;
    const auto outlet_static = outlet_evaluation.value;

    // 2. Independent inlet/outlet continuity.
    const IndependentFlowOperatingPoint flow{
        trial.inlet_mass_flow,
        trial.outlet_mass_flow,
        inlet_static.rho,
        outlet_static.rho
    };

    const auto continuity =
        evaluate_continuity_independent(
            geometry.inlet,
            geometry.exit,
            flow);

    // 3. Impeller prediction using the slip correlation.
    const ImpellerOperatingPoint operating_point{
        trial.rotational_speed_rpm,
        continuity.exit_meridional_velocity,
        trial.slip_factor,
        continuity.inlet_meridional_velocity,
        trial.inlet_prewhirl
    };

    const auto impeller =
        evaluate_impeller(
            geometry.impeller,
            operating_point);

    // 4. Actual trial velocities for thermodynamics.
    //
    // Importantly, outlet stagnation properties use
    // trial.outlet_whirl, NOT the slip prediction.
    const FlowStation inlet_station{
        inlet_static,
        {
            continuity.inlet_meridional_velocity,
            trial.inlet_prewhirl
        }
    };

    const FlowStation outlet_station{
        outlet_static,
        {
            continuity.exit_meridional_velocity,
            trial.outlet_whirl
        }
    };

    const auto inlet_stagnation =
        evaluate_stagnation_state(
            model,
            inlet_station);

    const auto outlet_stagnation =
        evaluate_stagnation_state(
            model,
            outlet_station);

    // 5. Independent conservation residuals.
    const double omega =
        2.0 * std::numbers::pi *
        trial.rotational_speed_rpm / 60.0;

    const ConservationGeometry conservation_geometry{
        geometry.impeller.inlet_diameter / 2.0,
        geometry.impeller.exit_diameter / 2.0
    };

    const ConservationState conservation_state{
        trial.inlet_mass_flow,
        trial.outlet_mass_flow,
        trial.inlet_prewhirl,
        trial.outlet_whirl,
        inlet_stagnation.stagnation_enthalpy,
        outlet_stagnation.stagnation_enthalpy,
        trial.shaft_torque,
        omega
    };

    const auto conservation =
        evaluate_conservation(
            conservation_geometry,
            conservation_state);

    // 6. Independent aerodynamic closure residual.
    const double slip_residual =
        trial.outlet_whirl -
        impeller.actual_exit_whirl;

    // Euler work calculated using the trial whirl,
    // rather than the slip-predicted whirl.
    const double euler_work =
        omega *
        (conservation_geometry.outlet_radius *
             trial.outlet_whirl -
         conservation_geometry.inlet_radius *
             trial.inlet_prewhirl);

    const double delta_h0 =
        outlet_stagnation.stagnation_enthalpy -
        inlet_stagnation.stagnation_enthalpy;

    if (!std::isfinite(slip_residual) ||
        !std::isfinite(euler_work) ||
        !std::isfinite(delta_h0)) {
        throw std::runtime_error(
            "Non-finite integrated impeller result");
    }

    return {
        inlet_static,
        outlet_static,
        continuity,
        impeller,
        inlet_stagnation,
        outlet_stagnation,
        conservation,
        impeller.actual_exit_whirl,
        trial.outlet_whirl,
        slip_residual,
        omega,
        euler_work,
        delta_h0
    };
}

} // namespace turbocompgpu
