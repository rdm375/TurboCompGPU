#include <turbocompgpu/impeller_system_loss.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-11)
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

thermomodel::PengRobinson methane_model()
{
    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(
            2200.0);

    return thermomodel::PengRobinson(caloric);
}

void test_ideal_gas_integration()
{
    const thermomodel::IdealGas model;

    const auto result =
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {20.0});

    const double expected_entropy_change =
        result.system.outlet_stagnation
            .stagnation_state.s
        -
        result.system.inlet_stagnation
            .stagnation_state.s;

    check_close(
        "Ideal-gas entropy change",
        result.loss.actual_entropy_increase,
        expected_entropy_change);

    check_close(
        "Ideal-gas loss residual",
        result.loss.entropy_residual,
        expected_entropy_change - 20.0);

    // Independent analytical ideal-gas identity.
    constexpr double cp = 1004.5;
    constexpr double R = 287.05;

    const auto& inlet =
        result.system.inlet_stagnation
            .stagnation_state;

    const auto& outlet =
        result.system.outlet_stagnation
            .stagnation_state;

    const double analytical_entropy_change =
        cp * std::log(outlet.T / inlet.T)
        -
        R * std::log(outlet.p / inlet.p);

    check_close(
        "Analytical stagnation entropy identity",
        result.loss.actual_entropy_increase,
        analytical_entropy_change);

    std::cout
        << "PASS: Integrated ideal-gas entropy loss\n";
}

void test_peng_robinson_integration()
{
    const auto model = methane_model();

    const auto result =
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {20.0});

    const auto& inlet =
        result.system.inlet_stagnation
            .stagnation_state;

    const auto& outlet =
        result.system.outlet_stagnation
            .stagnation_state;

    check_close(
        "PR entropy difference",
        result.loss.actual_entropy_increase,
        outlet.s - inlet.s);

    check_close(
        "PR entropy residual",
        result.loss.entropy_residual,
        outlet.s - inlet.s - 20.0);

    check_close(
        "PR conservation mass residual",
        result.system.conservation.mass,
        0.0);

    std::cout
        << "PASS: Integrated Peng-Robinson entropy loss\n";
}

void test_independent_loss_perturbation()
{
    const thermomodel::IdealGas model;

    const auto baseline =
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {10.0});

    const auto perturbed =
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {15.0});

    check_close(
        "Loss residual perturbation",
        perturbed.loss.entropy_residual -
            baseline.loss.entropy_residual,
        -5.0);

    check_close(
        "Mass residual invariance",
        perturbed.system.conservation.mass,
        baseline.system.conservation.mass);

    check_close(
        "Angular momentum residual invariance",
        perturbed.system.conservation.angular_momentum,
        baseline.system.conservation.angular_momentum);

    check_close(
        "Energy residual invariance",
        perturbed.system.conservation.energy,
        baseline.system.conservation.energy);

    check_close(
        "Slip residual invariance",
        perturbed.system.slip_residual,
        baseline.system.slip_residual);

    check_close(
        "Stagnation enthalpy invariance",
        perturbed.system.outlet_stagnation
            .stagnation_enthalpy,
        baseline.system.outlet_stagnation
            .stagnation_enthalpy);

    std::cout
        << "PASS: Independent loss perturbation\n";
}

void test_entropy_reference_invariance()
{
    const thermomodel::IdealGas model;

    const auto result =
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {20.0});

    auto inlet =
        result.system.inlet_stagnation
            .stagnation_state;

    auto outlet =
        result.system.outlet_stagnation
            .stagnation_state;

    constexpr double entropy_shift = 10000.0;

    inlet.s += entropy_shift;
    outlet.s += entropy_shift;

    const auto shifted =
        evaluate_entropy_loss(
            inlet,
            outlet,
            {20.0});

    check_close(
        "Entropy reference invariance",
        shifted.entropy_residual,
        result.loss.entropy_residual);

    std::cout
        << "PASS: Entropy reference invariance\n";
}

void test_invalid_loss()
{
    const thermomodel::IdealGas model;

    bool rejected = false;

    try {
        evaluate_impeller_system_with_loss(
            model,
            geometry(),
            trial(),
            {-1.0});
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    if (!rejected) {
        throw std::runtime_error(
            "Negative entropy generation accepted");
    }

    std::cout
        << "PASS: Invalid entropy-loss rejection\n";
}

} // namespace

int main()
{
    try {
        test_ideal_gas_integration();
        test_peng_robinson_integration();
        test_independent_loss_perturbation();
        test_entropy_reference_invariance();
        test_invalid_loss();

        std::cout
            << "\nPASS: All integrated entropy-loss tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
