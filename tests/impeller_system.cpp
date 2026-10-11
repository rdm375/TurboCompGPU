#include <turbocompgpu/impeller_system.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

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

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected << '\n';

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

ImpellerTrialState trial()
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

void test_integrated_evaluation()
{
    const thermomodel::IdealGas model;

    const auto g = geometry();
    const auto t = trial();

    const auto result =
        evaluate_impeller_system(model, g, t);

    check_close(
        "Inlet mass flow",
        result.continuity.inlet_mass_flow,
        t.inlet_mass_flow);

    check_close(
        "Outlet mass flow",
        result.continuity.exit_mass_flow,
        t.outlet_mass_flow);

    check_close(
        "Mass conservation",
        result.conservation.mass,
        0.0);

    check_close(
        "Slip residual",
        result.slip_residual,
        t.outlet_whirl -
            result.predicted_outlet_whirl);

    check_close(
        "Euler work",
        result.euler_specific_work,
        result.angular_speed *
        (g.impeller.exit_diameter / 2.0 *
             t.outlet_whirl -
         g.impeller.inlet_diameter / 2.0 *
             t.inlet_prewhirl));

    check_close(
        "Euler/slip work identity",
        result.euler_specific_work -
            result.impeller.euler_specific_work,
        result.impeller.tip_speed *
            result.slip_residual);

    check_close(
        "Inlet stagnation enthalpy",
        result.inlet_stagnation.stagnation_enthalpy,
        result.inlet_static.h +
            0.5 *
            (std::pow(
                 result.continuity.inlet_meridional_velocity,
                 2.0) +
             t.inlet_prewhirl * t.inlet_prewhirl));

    check_close(
        "Outlet stagnation enthalpy",
        result.outlet_stagnation.stagnation_enthalpy,
        result.outlet_static.h +
            0.5 *
            (std::pow(
                 result.continuity.exit_meridional_velocity,
                 2.0) +
             t.outlet_whirl * t.outlet_whirl));

    std::cout
        << "PASS: Integrated ideal-gas evaluation\n";
}

void test_independent_mass_residual()
{
    const thermomodel::IdealGas model;

    auto t = trial();
    t.outlet_mass_flow += 0.5;

    const auto result =
        evaluate_impeller_system(
            model, geometry(), t);

    check_close(
        "Independent mass residual",
        result.conservation.mass,
        0.5);

    std::cout
        << "PASS: Independent mass residual\n";
}

void test_independent_torque()
{
    const thermomodel::IdealGas model;

    const auto baseline =
        evaluate_impeller_system(
            model, geometry(), trial());

    auto perturbed_trial = trial();
    perturbed_trial.shaft_torque += 10.0;

    const auto perturbed =
        evaluate_impeller_system(
            model, geometry(), perturbed_trial);

    check_close(
        "Torque/angular momentum response",
        perturbed.conservation.angular_momentum -
            baseline.conservation.angular_momentum,
        10.0);

    check_close(
        "Torque/energy response",
        perturbed.conservation.energy -
            baseline.conservation.energy,
        -10.0 * baseline.angular_speed);

    check_close(
        "Torque does not change slip",
        perturbed.slip_residual,
        baseline.slip_residual);

    check_close(
        "Torque does not change stagnation enthalpy",
        perturbed.outlet_stagnation.stagnation_enthalpy,
        baseline.outlet_stagnation.stagnation_enthalpy);

    std::cout
        << "PASS: Independent torque perturbation\n";
}

void test_outlet_whirl()
{
    const thermomodel::IdealGas model;

    const auto baseline =
        evaluate_impeller_system(
            model, geometry(), trial());

    auto perturbed_trial = trial();
    perturbed_trial.outlet_whirl += 10.0;

    const auto perturbed =
        evaluate_impeller_system(
            model, geometry(), perturbed_trial);

    check_close(
        "Whirl/slip response",
        perturbed.slip_residual -
            baseline.slip_residual,
        10.0);

    check_close(
        "Whirl does not change slip prediction",
        perturbed.predicted_outlet_whirl,
        baseline.predicted_outlet_whirl);

    const double expected_delta_h0 =
        0.5 *
        (perturbed_trial.outlet_whirl *
             perturbed_trial.outlet_whirl -
         trial().outlet_whirl *
             trial().outlet_whirl);

    check_close(
        "Whirl/stagnation enthalpy response",
        perturbed.outlet_stagnation.stagnation_enthalpy -
            baseline.outlet_stagnation.stagnation_enthalpy,
        expected_delta_h0);

    std::cout
        << "PASS: Independent outlet-whirl perturbation\n";
}

} // namespace

int main()
{
    try {
        test_integrated_evaluation();
        test_independent_mass_residual();
        test_independent_torque();
        test_outlet_whirl();

        std::cout
            << "\nPASS: All integrated impeller tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
