#include <turbocompgpu/operating_point.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    using namespace turbocompgpu;

    const thermomodel::IdealGas model;

    const ImpellerSystemGeometry geometry {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearOperatingPoint operating {
        300.0, 500000.0, 5.0, 12000.0,
        0.90, 20.0, 20.0
    };

    const NonlinearResidualScales scales {
        50000.0, 100.0, 500000.0
    };

    const NonlinearUnknowns guess {
        340.0, std::log(700000.0 / scales.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto result = evaluate_operating_point(
        model, geometry, {operating, guess}, scales, options);

    if (!result.converged ||
        !result.independently_validated)
        throw std::runtime_error(
            "Operating point did not converge and validate");

    if (!result.basic_physics_valid)
        throw std::runtime_error(
            "Basic physical-consistency checks failed");

    if (result.operating_limits !=
        OperatingLimitStatus::NotEvaluated)
        throw std::runtime_error(
            "Operating limits were incorrectly classified");

    if (!(result.inlet_density > 0.0 &&
          result.outlet_density > 0.0))
        throw std::runtime_error(
            "Invalid thermodynamic densities");

    if (std::abs(result.outlet_temperature -
                 328.99653060808760) > 1.0e-6)
        throw std::runtime_error("Unexpected outlet temperature");

    if (std::abs(result.outlet_pressure -
                 644067.8864092113) > 1.0e-3)
        throw std::runtime_error("Unexpected outlet pressure");

    if (std::abs(result.pressure_ratio -
                 result.outlet_pressure / operating.inlet_pressure)
        > 1.0e-12)
        throw std::runtime_error("Incorrect pressure ratio");

    if (std::abs(result.shaft_power -
                 operating.mass_flow * result.specific_work)
        > 1.0e-8)
        throw std::runtime_error("Incorrect shaft power");

    const auto invalid = evaluate_operating_point(
        model, geometry, {operating, guess}, scales,
        options, -1.0);

    if (invalid.converged || invalid.independently_validated)
        throw std::runtime_error(
            "Invalid tolerance was accepted");

    std::cout << "PASS: operating-point API\n";
}
