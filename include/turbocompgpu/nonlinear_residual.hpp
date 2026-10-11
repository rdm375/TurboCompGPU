#pragma once

#include <turbocompgpu/impeller_system_loss.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Prescribed operating conditions for the reduced
// two-variable compressor nonlinear system.

struct NonlinearOperatingPoint {
    double inlet_temperature;
    double inlet_pressure;

    double mass_flow;
    double rotational_speed_rpm;

    double slip_factor;
    double inlet_prewhirl;

    double entropy_generation;
};

// Unknown outlet thermodynamic state.
//
// Log pressure is dimensionless:
//
//     log_pressure = log(p2 / pressure_scale)
//
// pressure_scale is a fixed numerical scale, not a
// thermodynamic reference convention.

struct NonlinearUnknowns {
    double outlet_temperature;
    double log_pressure;
};

struct NonlinearResidualScales {
    double specific_work;
    double specific_entropy;
    double pressure;
};

struct NonlinearResidualResult {
    NonlinearUnknowns unknowns;

    double outlet_pressure;

    double energy_residual;
    double entropy_residual;

    double normalized_energy;
    double normalized_entropy;

    double residual_norm;

    double outlet_meridional_velocity;
    double outlet_whirl;
    double euler_specific_work;
    double shaft_torque;

    ImpellerSystemWithLossResult full_system;
};

inline NonlinearResidualResult evaluate_nonlinear_residual(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const NonlinearOperatingPoint& operating,
    const NonlinearUnknowns& unknowns,
    const NonlinearResidualScales& scales)
{
    if (!(std::isfinite(scales.specific_work) &&
          scales.specific_work > 0.0 &&
          std::isfinite(scales.specific_entropy) &&
          scales.specific_entropy > 0.0 &&
          std::isfinite(scales.pressure) &&
          scales.pressure > 0.0)) {
        throw std::invalid_argument(
            "Nonlinear residual scales must be positive and finite");
    }

    if (!(std::isfinite(unknowns.outlet_temperature) &&
          unknowns.outlet_temperature > 0.0 &&
          std::isfinite(unknowns.log_pressure))) {
        throw std::invalid_argument(
            "Invalid nonlinear unknowns");
    }

    const double outlet_pressure =
        scales.pressure * std::exp(unknowns.log_pressure);

    if (!(std::isfinite(outlet_pressure) &&
          outlet_pressure > 0.0)) {
        throw std::invalid_argument(
            "Nonlinear outlet pressure is invalid");
    }

    // Evaluate the outlet state once to determine
    // continuity velocity and slip-predicted whirl.

    const auto outlet_state_result =
        model.evaluate_Tp(
            unknowns.outlet_temperature,
            outlet_pressure);

    if (!outlet_state_result) {
        throw std::runtime_error(
            "Outlet thermodynamic evaluation failed: " +
            outlet_state_result.message);
    }

    const auto& outlet = outlet_state_result.value;

    if (!(std::isfinite(outlet.rho) &&
          outlet.rho > 0.0)) {
        throw std::runtime_error(
            "Invalid outlet density");
    }

    const double pi =
        3.141592653589793238462643383279502884;

    const double area =
        pi *
        geometry.exit.diameter *
        geometry.exit.width *
        (1.0 - geometry.exit.blockage);

    if (!(std::isfinite(area) && area > 0.0)) {
        throw std::invalid_argument(
            "Invalid outlet effective flow area");
    }

    const double Vm2 =
        operating.mass_flow / (outlet.rho * area);

    const double omega =
        2.0 * pi * operating.rotational_speed_rpm / 60.0;

    const double U2 =
        omega * geometry.impeller.exit_diameter / 2.0;

    const double U1 =
        omega * geometry.impeller.inlet_diameter / 2.0;

    const double beta =
        geometry.impeller.exit_blade_angle_deg *
        pi / 180.0;

    const double whirl =
        operating.slip_factor *
        (U2 - Vm2 / std::tan(beta));

    const double euler_work =
        U2 * whirl - U1 * operating.inlet_prewhirl;

    const double torque =
        operating.mass_flow *
        (geometry.impeller.exit_diameter * whirl / 2.0 -
         geometry.impeller.inlet_diameter *
             operating.inlet_prewhirl / 2.0);

    // The full system remains our independent conservation
    // and thermodynamic validation machinery.

    const ImpellerTrialState trial {
        operating.inlet_temperature,
        operating.inlet_pressure,
        unknowns.outlet_temperature,
        outlet_pressure,
        operating.mass_flow,
        operating.mass_flow,
        operating.rotational_speed_rpm,
        operating.slip_factor,
        operating.inlet_prewhirl,
        whirl,
        torque
    };

    const auto full =
        evaluate_impeller_system_with_loss(
            model,
            geometry,
            trial,
            {operating.entropy_generation});

    const double h01 =
        full.system.inlet_stagnation.stagnation_enthalpy;

    const double h02 =
        full.system.outlet_stagnation.stagnation_enthalpy;

    const double energy_residual =
        (h02 - h01) - euler_work;

    const double entropy_residual =
        full.loss.entropy_residual;

    const double normalized_energy =
        energy_residual / scales.specific_work;

    const double normalized_entropy =
        entropy_residual / scales.specific_entropy;

    const double residual_norm =
        std::hypot(normalized_energy, normalized_entropy);

    if (!std::isfinite(residual_norm)) {
        throw std::runtime_error(
            "Non-finite nonlinear residual");
    }

    return {
        unknowns,
        outlet_pressure,
        energy_residual,
        entropy_residual,
        normalized_energy,
        normalized_entropy,
        residual_norm,
        Vm2,
        whirl,
        euler_work,
        torque,
        full
    };
}

} // namespace turbocompgpu
