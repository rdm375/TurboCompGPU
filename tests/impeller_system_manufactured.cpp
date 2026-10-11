#include <turbocompgpu/impeller_system.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-10)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    const double error =
        std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) ||
        error > tolerance) {
        std::cerr
            << std::setprecision(17)
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected
            << " relative_error=" << error
            << '\n';

        throw std::runtime_error(name);
    }
}

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

ImpellerTrialState baseline_trial()
{
    return {
        300.0,
        500000.0,

        340.0,
        700000.0,

        5.0,
        5.0,

        12000.0,
        0.90,

        20.0,
        180.0,

        200.0
    };
}

thermomodel::PengRobinson methane_model()
{
    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(
            2200.0);

    return thermomodel::PengRobinson(caloric);
}

void test_peng_robinson_integration()
{
    const auto model = methane_model();

    const auto g = geometry();
    const auto trial = baseline_trial();

    const auto result =
        evaluate_impeller_system(
            model, g, trial);

    if (!(result.inlet_static.rho > 0.0) ||
        !(result.outlet_static.rho > 0.0)) {
        throw std::runtime_error(
            "Nonpositive Peng-Robinson density");
    }

    check_close(
        "PR inlet mass flow",
        result.continuity.inlet_mass_flow,
        trial.inlet_mass_flow);

    check_close(
        "PR outlet mass flow",
        result.continuity.exit_mass_flow,
        trial.outlet_mass_flow);

    check_close(
        "PR mass residual",
        result.conservation.mass,
        0.0);

    check_close(
        "PR slip residual",
        result.slip_residual,
        trial.outlet_whirl -
            result.predicted_outlet_whirl);

    check_close(
        "PR Euler/slip identity",
        result.euler_specific_work -
            result.impeller.euler_specific_work,
        result.impeller.tip_speed *
            result.slip_residual);

    check_close(
        "PR inlet stagnation enthalpy",
        result.inlet_stagnation.stagnation_enthalpy,
        result.inlet_static.h +
            0.5 *
            (std::pow(
                 result.continuity.inlet_meridional_velocity,
                 2.0) +
             trial.inlet_prewhirl *
                 trial.inlet_prewhirl));

    check_close(
        "PR outlet stagnation enthalpy",
        result.outlet_stagnation.stagnation_enthalpy,
        result.outlet_static.h +
            0.5 *
            (std::pow(
                 result.continuity.exit_meridional_velocity,
                 2.0) +
             trial.outlet_whirl *
                 trial.outlet_whirl));

    std::cout
        << "PASS: Peng-Robinson integrated evaluation\n";
}

// Evaluate the scalar energy/slip consistency equation
// without invoking the integrated impeller-system evaluator.
//
// This provides an independent construction route for
// the manufactured conservative operating point.

double manufactured_energy_defect(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& g,
    const ImpellerTrialState& trial,
    double outlet_temperature)
{
    const auto inlet_result =
        model.evaluate_Tp(
            trial.inlet_temperature,
            trial.inlet_pressure);

    const auto outlet_result =
        model.evaluate_Tp(
            outlet_temperature,
            trial.outlet_pressure);

    if (!inlet_result || !outlet_result) {
        throw std::runtime_error(
            "Manufactured-state thermodynamic evaluation failed");
    }

    const auto& inlet = inlet_result.value;
    const auto& outlet = outlet_result.value;

    const double A1 =
        effective_flow_area(
            inlet_annulus_area(g.inlet),
            g.inlet.blockage);

    const double A2 =
        effective_flow_area(
            impeller_exit_area(g.exit),
            g.exit.blockage);

    const double Vm1 =
        trial.inlet_mass_flow /
        (inlet.rho * A1);

    const double Vm2 =
        trial.outlet_mass_flow /
        (outlet.rho * A2);

    const double omega =
        2.0 * std::acos(-1.0) *
        trial.rotational_speed_rpm / 60.0;

    const double r1 =
        g.impeller.inlet_diameter / 2.0;

    const double r2 =
        g.impeller.exit_diameter / 2.0;

    const double U1 = omega * r1;
    const double U2 = omega * r2;

    const double beta =
        g.impeller.exit_blade_angle_deg *
        std::acos(-1.0) / 180.0;

    const double predicted_whirl =
        trial.slip_factor *
        (U2 - Vm2 / std::tan(beta));

    const double h01 =
        inlet.h +
        0.5 *
        (Vm1 * Vm1 +
         trial.inlet_prewhirl *
             trial.inlet_prewhirl);

    const double h02 =
        outlet.h +
        0.5 *
        (Vm2 * Vm2 +
         predicted_whirl * predicted_whirl);

    const double euler_work =
        U2 * predicted_whirl -
        U1 * trial.inlet_prewhirl;

    return h02 - h01 - euler_work;
}

double solve_manufactured_temperature(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& g,
    const ImpellerTrialState& trial)
{
    // Bracket the outlet-temperature root.
    // Both endpoints are inside the configured model domain.

    double lower = 200.0;
    double upper = 600.0;

    double f_lower =
        manufactured_energy_defect(
            model, g, trial, lower);

    double f_upper =
        manufactured_energy_defect(
            model, g, trial, upper);

    if (!std::isfinite(f_lower) ||
        !std::isfinite(f_upper) ||
        std::signbit(f_lower) ==
            std::signbit(f_upper)) {
        throw std::runtime_error(
            "Manufactured-state root not bracketed");
    }

    for (int iteration = 0;
         iteration < 100;
         ++iteration) {

        const double midpoint =
            lower + 0.5 * (upper - lower);

        const double f_midpoint =
            manufactured_energy_defect(
                model, g, trial, midpoint);

        if (!std::isfinite(f_midpoint)) {
            throw std::runtime_error(
                "Non-finite manufactured-state defect");
        }

        if (std::signbit(f_midpoint) ==
            std::signbit(f_lower)) {
            lower = midpoint;
            f_lower = f_midpoint;
        }
        else {
            upper = midpoint;
            f_upper = f_midpoint;
        }

        if (upper - lower <=
            1.0e-12 *
            std::max(1.0, std::abs(midpoint))) {
            return lower +
                0.5 * (upper - lower);
        }
    }

    throw std::runtime_error(
        "Manufactured-state temperature did not converge");
}

void test_manufactured_state(
    const thermomodel::Model& model,
    const char* label)
{
    const auto g = geometry();

    auto trial = baseline_trial();

    // Equal inlet/outlet mass flows are necessary
    // for a steady conservative manufactured state.
    trial.outlet_mass_flow =
        trial.inlet_mass_flow;

    trial.outlet_temperature =
        solve_manufactured_temperature(
            model, g, trial);

    // Determine slip-predicted whirl from the
    // independently evaluated outlet density.
    const auto outlet_result =
        model.evaluate_Tp(
            trial.outlet_temperature,
            trial.outlet_pressure);

    if (!outlet_result) {
        throw std::runtime_error(
            "Manufactured outlet evaluation failed");
    }

    const double A2 =
        effective_flow_area(
            impeller_exit_area(g.exit),
            g.exit.blockage);

    const double Vm2 =
        trial.outlet_mass_flow /
        (outlet_result.value.rho * A2);

    const double omega =
        2.0 * std::acos(-1.0) *
        trial.rotational_speed_rpm / 60.0;

    const double r1 =
        g.impeller.inlet_diameter / 2.0;

    const double r2 =
        g.impeller.exit_diameter / 2.0;

    const double U2 = omega * r2;

    const double beta =
        g.impeller.exit_blade_angle_deg *
        std::acos(-1.0) / 180.0;

    trial.outlet_whirl =
        trial.slip_factor *
        (U2 - Vm2 / std::tan(beta));

    trial.shaft_torque =
        trial.inlet_mass_flow *
        (r2 * trial.outlet_whirl -
         r1 * trial.inlet_prewhirl);

    const auto result =
        evaluate_impeller_system(
            model, g, trial);

    const double mass_scale =
        trial.inlet_mass_flow;

    const double torque_scale =
        std::max(
            1.0,
            std::abs(trial.shaft_torque));

    const double power_scale =
        std::max(
            1.0,
            std::abs(
                trial.shaft_torque * omega));

    check_close(
        "Manufactured mass residual",
        result.conservation.mass /
            mass_scale,
        0.0,
        1.0e-12);

    check_close(
        "Manufactured angular-momentum residual",
        result.conservation.angular_momentum /
            torque_scale,
        0.0,
        1.0e-12);

    check_close(
        "Manufactured energy residual",
        result.conservation.energy /
            power_scale,
        0.0,
        1.0e-10);

    check_close(
        "Manufactured slip residual",
        result.slip_residual /
            std::max(
                1.0,
                std::abs(trial.outlet_whirl)),
        0.0,
        1.0e-12);

    std::cout
        << std::setprecision(17)
        << label
        << " manufactured state:\n"
        << "  T2             = "
        << trial.outlet_temperature << " K\n"
        << "  Vtheta2        = "
        << trial.outlet_whirl << " m/s\n"
        << "  shaft torque   = "
        << trial.shaft_torque << " Nm\n"
        << "  mass residual  = "
        << result.conservation.mass << " kg/s\n"
        << "  torque residual= "
        << result.conservation.angular_momentum
        << " Nm\n"
        << "  energy residual= "
        << result.conservation.energy << " W\n"
        << "  slip residual  = "
        << result.slip_residual << " m/s\n";

    std::cout
        << "PASS: " << label
        << " manufactured conservation state\n";
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal_gas;
        const auto methane = methane_model();

        test_peng_robinson_integration();

        test_manufactured_state(
            ideal_gas,
            "Ideal gas");

        test_manufactured_state(
            methane,
            "Peng-Robinson methane");

        std::cout
            << "\nPASS: All manufactured-state tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
