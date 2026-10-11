#include <turbocompgpu/station.hpp>

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

    if (!std::isfinite(actual) || error > tolerance) {
        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected
            << " relative_error=" << error
            << '\n';

        throw std::runtime_error(name);
    }
}

template<class Function>
void expect_failure(const char* name, Function&& function)
{
    try {
        function();
    }
    catch (const std::exception&) {
        std::cout << "PASS: " << name << '\n';
        return;
    }

    throw std::runtime_error(name);
}

void test_ideal_gas()
{
    using namespace turbocompgpu;

    thermomodel::IdealGas model;

    const auto static_state =
        evaluate(model, 300.0, 5.0e6);

    const FlowStation station{
        static_state,
        {120.0, 160.0}
    };

    const auto result =
        evaluate_stagnation_state(model, station);

    const auto& gas = model.parameters();

    const double V2 =
        120.0 * 120.0 + 160.0 * 160.0;

    const double T0 =
        300.0 + V2 / (2.0 * gas.cp);

    const double p0 =
        5.0e6 * std::pow(
            T0 / 300.0,
            gas.cp / gas.R);

    check_close(
        "Ideal gas stagnation temperature",
        result.stagnation_state.T,
        T0);

    check_close(
        "Ideal gas stagnation pressure",
        result.stagnation_state.p,
        p0);

    check_close(
        "Ideal gas stagnation enthalpy",
        result.stagnation_enthalpy,
        static_state.h + 0.5 * V2);

    check_close(
        "Ideal gas entropy conservation",
        result.stagnation_state.s,
        static_state.s);

    // A rotation of the absolute velocity vector must
    // leave the stagnation state unchanged.
    const FlowStation rotated{
        static_state,
        {160.0, -120.0}
    };

    const auto rotated_result =
        evaluate_stagnation_state(model, rotated);

    check_close(
        "Velocity direction invariance T0",
        rotated_result.stagnation_state.T,
        result.stagnation_state.T);

    check_close(
        "Velocity direction invariance p0",
        rotated_result.stagnation_state.p,
        result.stagnation_state.p);

    std::cout
        << "Ideal gas T0 = "
        << result.stagnation_state.T << " K\n"
        << "Ideal gas p0 = "
        << result.stagnation_state.p << " Pa\n"
        << "Pressure iterations = "
        << result.pressure_iterations << '\n'
        << "Thermodynamic evaluations = "
        << result.temperature_evaluations << '\n';

    std::cout
        << "PASS: Ideal gas stagnation validation\n";
}

void test_zero_velocity()
{
    using namespace turbocompgpu;

    thermomodel::IdealGas model;

    const auto static_state =
        evaluate(model, 300.0, 5.0e6);

    const auto result =
        evaluate_stagnation_state(
            model,
            {static_state, {0.0, 0.0}});

    check_close(
        "Zero velocity T0",
        result.stagnation_state.T,
        static_state.T);

    check_close(
        "Zero velocity p0",
        result.stagnation_state.p,
        static_state.p);

    if (result.pressure_evaluations != 0) {
        throw std::runtime_error(
            "Zero velocity unexpectedly required pressure inversion");
    }

    std::cout
        << "PASS: Zero velocity stagnation state\n";
}

void test_peng_robinson()
{
    using namespace turbocompgpu;

    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(2200.0);

    thermomodel::PRParameters params;
    params.root_policy =
        thermomodel::PRRootPolicy::Vapor;

    thermomodel::PengRobinson model(
        caloric, params);

    const auto static_state =
        evaluate(model, 300.0, 5.0e6);

    const FlowStation station{
        static_state,
        {100.0, 150.0}
    };

    const auto result =
        evaluate_stagnation_state(model, station);

    const double kinetic_energy =
        0.5 * (100.0 * 100.0 + 150.0 * 150.0);

    check_close(
        "PR stagnation energy",
        result.stagnation_state.h - static_state.h,
        kinetic_energy);

    check_close(
        "PR stagnation entropy",
        result.stagnation_state.s,
        static_state.s);

    if (!(result.stagnation_state.T > static_state.T) ||
        !(result.stagnation_state.p > static_state.p)) {
        throw std::runtime_error(
            "Unexpected PR stagnation state ordering");
    }

    std::cout
        << "PR T0 = "
        << result.stagnation_state.T << " K\n"
        << "PR p0 = "
        << result.stagnation_state.p << " Pa\n"
        << "PR enthalpy residual = "
        << result.enthalpy_residual << " J/kg\n"
        << "PR entropy residual = "
        << result.entropy_residual << " J/(kg K)\n"
        << "Pressure iterations = "
        << result.pressure_iterations << '\n'
        << "Thermodynamic evaluations = "
        << result.temperature_evaluations << '\n';

    std::cout
        << "PASS: Peng-Robinson stagnation validation\n";
}

void test_invalid_velocity()
{
    using namespace turbocompgpu;

    thermomodel::IdealGas model;

    const auto state =
        evaluate(model, 300.0, 5.0e6);

    const double nan =
        std::numeric_limits<double>::quiet_NaN();

    const double inf =
        std::numeric_limits<double>::infinity();

    expect_failure("NaN velocity", [&] {
        evaluate_stagnation_state(
            model, {state, {nan, 0.0}});
    });

    expect_failure("Infinite velocity", [&] {
        evaluate_stagnation_state(
            model, {state, {0.0, inf}});
    });
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(17);

        test_ideal_gas();
        test_zero_velocity();
        test_peng_robinson();
        test_invalid_velocity();

        std::cout
            << "PASS: All station thermodynamics tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
