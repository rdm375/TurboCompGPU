#include <turbocompgpu/choking_mass_flux.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void check_close(
    double actual,
    double expected,
    double tolerance,
    const char* message)
{
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            tolerance * std::max(1.0, std::abs(expected))) {
        throw std::runtime_error(message);
    }
}

void test_ideal_gas_choking()
{
    thermomodel::IdealGas model;

    const auto& params = model.parameters();

    const double gamma =
        params.cp / (params.cp - params.R);

    constexpr double T0 = 300.0;
    constexpr double p0 = 500000.0;

    const double expected_pressure =
        p0 * std::pow(
            2.0 / (gamma + 1.0),
            gamma / (gamma - 1.0));

    const double expected_flux =
        p0 / std::sqrt(T0) *
        std::sqrt(gamma / params.R) *
        std::pow(
            2.0 / (gamma + 1.0),
            (gamma + 1.0) /
            (2.0 * (gamma - 1.0)));

    const auto result =
        turbocompgpu::evaluate_choking_mass_flux(
            model, model, T0, p0);

    check_close(
        result.pressure,
        expected_pressure,
        1.0e-8,
        "Ideal-gas critical pressure mismatch");

    check_close(
        result.critical_mass_flux,
        expected_flux,
        1.0e-8,
        "Ideal-gas critical mass flux mismatch");

    check_close(
        result.entropy_residual,
        0.0,
        1.0e-8,
        "Ideal-gas entropy residual");

    std::cout
        << "PASS: ideal-gas analytical choking oracle\n";
}

void test_peng_robinson_choking()
{
    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(2200.0);

    thermomodel::PRParameters parameters;
    parameters.root_policy =
        thermomodel::PRRootPolicy::Vapor;

    thermomodel::PengRobinson model(
        caloric, parameters);

    constexpr double T0 = 300.0;
    constexpr double p0 = 500000.0;

    const auto result =
        turbocompgpu::evaluate_choking_mass_flux(
            model, model, T0, p0);

    if (!(result.critical_mass_flux > 0.0 &&
          result.sound_speed > 0.0 &&
          result.pressure < p0)) {
        throw std::runtime_error(
            "Invalid Peng-Robinson choking state");
    }

    check_close(
        result.entropy_residual,
        0.0,
        1.0e-8,
        "Peng-Robinson entropy residual");

    std::cout
        << std::setprecision(12)
        << "PR critical pressure: "
        << result.pressure << " Pa\n"
        << "PR critical temperature: "
        << result.temperature << " K\n"
        << "PR critical mass flux: "
        << result.critical_mass_flux << " kg/(m2 s)\n"
        << "PR sonic speed: "
        << result.sound_speed << " m/s\n";

    std::cout
        << "PASS: Peng-Robinson choking calculation\n";
}

} // namespace

int main()
{
    try {
        test_ideal_gas_choking();
        test_peng_robinson_choking();

        std::cout
            << "PASS: All choking mass-flux tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
