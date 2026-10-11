#include <turbocompgpu/nonlinear_residual.hpp>
#include <turbocompgpu/reference_shift_model.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
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

NonlinearOperatingPoint operating()
{
    return {
        300.0,
        500000.0,
        5.0,
        12000.0,
        0.90,
        20.0,
        20.0
    };
}

NonlinearResidualScales scales()
{
    return {
        50000.0,
        100.0,
        500000.0
    };
}

void test_manufactured_solution()
{
    const thermomodel::IdealGas model;

    const NonlinearUnknowns unknowns {
        328.9965306080876,
        std::log(644067.88640921132 / scales().pressure)
    };

    const auto result =
        evaluate_nonlinear_residual(
            model,
            geometry(),
            operating(),
            unknowns,
            scales());

    require(
        std::abs(result.normalized_energy) < 1.0e-10,
        "Manufactured energy residual");

    require(
        std::abs(result.normalized_entropy) < 1.0e-10,
        "Manufactured entropy residual");

    require(
        std::abs(result.full_system.system.conservation.mass)
            < 1.0e-12,
        "Manufactured mass residual");

    require(
        std::abs(result.full_system.system.conservation.angular_momentum)
            < 1.0e-10,
        "Manufactured angular momentum residual");

    require(
        std::abs(result.full_system.system.slip_residual)
            < 1.0e-10,
        "Manufactured slip residual");

    std::cout
        << "PASS: Manufactured nonlinear residual\n";
}

void test_independent_perturbations()
{
    const thermomodel::IdealGas model;

    const NonlinearUnknowns base {
        328.9965306080876,
        std::log(644067.88640921132 / scales().pressure)
    };

    const auto baseline =
        evaluate_nonlinear_residual(
            model, geometry(), operating(), base, scales());

    auto temperature = base;
    temperature.outlet_temperature += 1.0;

    const auto changed_temperature =
        evaluate_nonlinear_residual(
            model, geometry(), operating(),
            temperature, scales());

    auto pressure = base;
    pressure.log_pressure += 0.01;

    const auto changed_pressure =
        evaluate_nonlinear_residual(
            model, geometry(), operating(),
            pressure, scales());

    require(
        std::abs(changed_temperature.normalized_energy -
                 baseline.normalized_energy) > 1.0e-4,
        "Temperature perturbation did not affect energy");

    require(
        std::abs(changed_pressure.normalized_entropy -
                 baseline.normalized_entropy) > 1.0e-4,
        "Pressure perturbation did not affect entropy");

    std::cout
        << "PASS: Independent nonlinear perturbations\n";
}

void test_reference_invariance()
{
    const thermomodel::IdealGas model;

    const ReferenceShiftModel shifted(
        model, 1.0e6, 1.0e4);

    const NonlinearUnknowns unknowns {
        328.9965306080876,
        std::log(644067.88640921132 / scales().pressure)
    };

    const auto original =
        evaluate_nonlinear_residual(
            model, geometry(), operating(),
            unknowns, scales());

    const auto modified =
        evaluate_nonlinear_residual(
            shifted, geometry(), operating(),
            unknowns, scales());

    require(
        std::abs(original.normalized_energy -
                 modified.normalized_energy) < 1.0e-10,
        "Energy residual depends on enthalpy reference");

    require(
        std::abs(original.normalized_entropy -
                 modified.normalized_entropy) < 1.0e-10,
        "Entropy residual depends on entropy reference");

    require(
        std::abs(original.outlet_pressure -
                 modified.outlet_pressure) < 1.0e-9,
        "Outlet pressure depends on reference convention");

    require(
        std::abs(original.outlet_whirl -
                 modified.outlet_whirl) < 1.0e-12,
        "Outlet whirl depends on reference convention");

    std::cout
        << "PASS: Nonlinear reference invariance\n";
}

} // namespace

int main()
{
    try {
        test_manufactured_solution();
        test_independent_perturbations();
        test_reference_invariance();

        std::cout
            << "PASS: All nonlinear residual tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
