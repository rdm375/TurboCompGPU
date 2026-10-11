#pragma once

#include <turbocompgpu/nonlinear_residual.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

struct FastResidualContext {
    const thermomodel::Model* model = nullptr;

    ImpellerSystemGeometry geometry;
    NonlinearOperatingPoint operating;
    NonlinearResidualScales scales;

    thermomodel::State inlet_state;

    double inlet_stagnation_enthalpy = 0.0;
    double inlet_entropy = 0.0;

    double inlet_meridional_velocity = 0.0;
    double inlet_blade_speed = 0.0;
    double outlet_blade_speed = 0.0;
    double outlet_effective_area = 0.0;

    double angular_speed = 0.0;
    double cotangent_exit_angle = 0.0;
};

inline FastResidualContext make_fast_residual_context(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const NonlinearOperatingPoint& operating,
    const NonlinearResidualScales& scales)
{
    constexpr double pi =
        3.141592653589793238462643383279502884;

    if (!(std::isfinite(operating.mass_flow) &&
          operating.mass_flow > 0.0 &&
          std::isfinite(operating.rotational_speed_rpm) &&
          operating.rotational_speed_rpm > 0.0 &&
          std::isfinite(operating.slip_factor) &&
          operating.slip_factor > 0.0 &&
          operating.slip_factor <= 1.0 &&
          std::isfinite(operating.inlet_prewhirl) &&
          std::isfinite(operating.entropy_generation) &&
          operating.entropy_generation >= 0.0)) {
        throw std::invalid_argument(
            "Invalid fast residual operating point");
    }

    if (!(std::isfinite(scales.specific_work) &&
          scales.specific_work > 0.0 &&
          std::isfinite(scales.specific_entropy) &&
          scales.specific_entropy > 0.0 &&
          std::isfinite(scales.pressure) &&
          scales.pressure > 0.0)) {
        throw std::invalid_argument(
            "Invalid fast residual scales");
    }

    const auto& inlet = geometry.inlet;
    const auto& exit = geometry.exit;
    const auto& impeller = geometry.impeller;

    if (!(std::isfinite(inlet.hub_diameter) &&
          std::isfinite(inlet.shroud_diameter) &&
          inlet.hub_diameter >= 0.0 &&
          inlet.shroud_diameter > inlet.hub_diameter &&
          std::isfinite(inlet.blockage) &&
          inlet.blockage >= 0.0 &&
          inlet.blockage < 1.0 &&
          std::isfinite(exit.diameter) &&
          exit.diameter > 0.0 &&
          std::isfinite(exit.width) &&
          exit.width > 0.0 &&
          std::isfinite(exit.blockage) &&
          exit.blockage >= 0.0 &&
          exit.blockage < 1.0 &&
          std::isfinite(impeller.exit_diameter) &&
          impeller.exit_diameter == exit.diameter &&
          std::isfinite(impeller.inlet_diameter) &&
          impeller.inlet_diameter > 0.0 &&
          impeller.inlet_diameter < impeller.exit_diameter &&
          std::isfinite(impeller.exit_blade_angle_deg) &&
          impeller.exit_blade_angle_deg > 0.0 &&
          impeller.exit_blade_angle_deg <= 90.0)) {
        throw std::invalid_argument(
            "Invalid fast residual geometry");
    }

    const auto inlet_result = model.evaluate_Tp(
        operating.inlet_temperature,
        operating.inlet_pressure);

    if (!inlet_result) {
        throw std::runtime_error(
            "Fast residual inlet evaluation failed: " +
            inlet_result.message);
    }

    const auto& state = inlet_result.value;

    if (!(std::isfinite(state.rho) && state.rho > 0.0 &&
          std::isfinite(state.h) &&
          std::isfinite(state.s))) {
        throw std::runtime_error(
            "Invalid inlet thermodynamic properties");
    }

    const double inlet_area =
        (pi / 4.0) *
        (inlet.shroud_diameter * inlet.shroud_diameter -
         inlet.hub_diameter * inlet.hub_diameter) *
        (1.0 - inlet.blockage);

    const double outlet_area =
        pi * exit.diameter * exit.width *
        (1.0 - exit.blockage);

    if (!(std::isfinite(inlet_area) && inlet_area > 0.0 &&
          std::isfinite(outlet_area) && outlet_area > 0.0)) {
        throw std::invalid_argument(
            "Invalid effective flow area");
    }

    const double vm1 =
        operating.mass_flow / (state.rho * inlet_area);

    const double h01 =
        state.h +
        0.5 * (vm1 * vm1 +
               operating.inlet_prewhirl *
                   operating.inlet_prewhirl);

    const double omega =
        2.0 * pi * operating.rotational_speed_rpm / 60.0;

    const double beta =
        impeller.exit_blade_angle_deg * pi / 180.0;

    FastResidualContext context;

    context.model = &model;
    context.geometry = geometry;
    context.operating = operating;
    context.scales = scales;
    context.inlet_state = state;

    context.inlet_stagnation_enthalpy = h01;
    context.inlet_entropy = state.s;
    context.inlet_meridional_velocity = vm1;

    context.inlet_blade_speed =
        omega * impeller.inlet_diameter / 2.0;

    context.outlet_blade_speed =
        omega * impeller.exit_diameter / 2.0;

    context.outlet_effective_area = outlet_area;
    context.angular_speed = omega;

    context.cotangent_exit_angle =
        std::cos(beta) / std::sin(beta);

    return context;
}

inline NonlinearResidualResult evaluate_fast_nonlinear_residual(
    const FastResidualContext& context,
    const NonlinearUnknowns& unknowns)
{
    if (context.model == nullptr) {
        throw std::invalid_argument(
            "Fast residual context has no model");
    }

    if (!(std::isfinite(unknowns.outlet_temperature) &&
          unknowns.outlet_temperature > 0.0 &&
          std::isfinite(unknowns.log_pressure))) {
        throw std::invalid_argument(
            "Invalid fast residual unknowns");
    }

    const double pressure =
        context.scales.pressure *
        std::exp(unknowns.log_pressure);

    if (!(std::isfinite(pressure) && pressure > 0.0)) {
        throw std::invalid_argument(
            "Invalid fast residual outlet pressure");
    }

    const auto outlet_result =
        context.model->evaluate_Tp(
            unknowns.outlet_temperature,
            pressure);

    if (!outlet_result) {
        throw std::runtime_error(
            "Fast residual outlet evaluation failed: " +
            outlet_result.message);
    }

    const auto& outlet = outlet_result.value;

    if (!(std::isfinite(outlet.rho) &&
          outlet.rho > 0.0 &&
          std::isfinite(outlet.h) &&
          std::isfinite(outlet.s))) {
        throw std::runtime_error(
            "Invalid outlet thermodynamic properties");
    }

    const double vm2 =
        context.operating.mass_flow /
        (outlet.rho * context.outlet_effective_area);

    const double whirl =
        context.operating.slip_factor *
        (context.outlet_blade_speed -
         vm2 * context.cotangent_exit_angle);

    const double euler_work =
        context.outlet_blade_speed * whirl -
        context.inlet_blade_speed *
            context.operating.inlet_prewhirl;

    const double h02 =
        outlet.h +
        0.5 * (vm2 * vm2 + whirl * whirl);

    const double energy_residual =
        (h02 - context.inlet_stagnation_enthalpy) -
        euler_work;

    const double entropy_residual =
        (outlet.s - context.inlet_entropy) -
        context.operating.entropy_generation;

    const double normalized_energy =
        energy_residual / context.scales.specific_work;

    const double normalized_entropy =
        entropy_residual / context.scales.specific_entropy;

    const double norm =
        std::hypot(normalized_energy, normalized_entropy);

    const double torque =
        context.operating.mass_flow *
        (context.geometry.exit.diameter * whirl / 2.0 -
         context.geometry.impeller.inlet_diameter *
             context.operating.inlet_prewhirl / 2.0);

    if (!(std::isfinite(norm) &&
          std::isfinite(vm2) &&
          std::isfinite(whirl) &&
          std::isfinite(euler_work) &&
          std::isfinite(torque))) {
        throw std::runtime_error(
            "Non-finite fast residual result");
    }

    // The final member is intentionally default-initialized.
    //
    // Fast residual evaluation does NOT reconstruct or
    // validate the complete five-equation system.
    //
    // Consumers must explicitly call the original full
    // evaluator when independent validation is required.

    return {
        unknowns,
        pressure,
        energy_residual,
        entropy_residual,
        normalized_energy,
        normalized_entropy,
        norm,
        vm2,
        whirl,
        euler_work,
        torque,
        {}
    };
}

inline ImpellerSystemWithLossResult validate_fast_solution(
    const FastResidualContext& context,
    const NonlinearResidualResult& solution)
{
    if (context.model == nullptr) {
        throw std::invalid_argument(
            "Fast residual context has no model");
    }

    const ImpellerTrialState trial {
        context.operating.inlet_temperature,
        context.operating.inlet_pressure,
        solution.unknowns.outlet_temperature,
        solution.outlet_pressure,
        context.operating.mass_flow,
        context.operating.mass_flow,
        context.operating.rotational_speed_rpm,
        context.operating.slip_factor,
        context.operating.inlet_prewhirl,
        solution.outlet_whirl,
        solution.shaft_torque
    };

    return evaluate_impeller_system_with_loss(
        *context.model,
        context.geometry,
        trial,
        {context.operating.entropy_generation});
}

} // namespace turbocompgpu
