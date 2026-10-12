#include <turbocompgpu/compressor_map.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-9)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance * scale) {
        std::cerr << "FAIL: " << name
                  << " actual=" << actual
                  << " expected=" << expected << '\n';

        throw std::runtime_error(name);
    }
}

template<class Model>
void test_model(const Model& model, const char* name)
{
    const ImpellerSystemGeometry geometry {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearOperatingPoint operating {
        300.0, 500000.0,
        5.0, 12000.0,
        0.90, 20.0, 20.0
    };

    const NonlinearResidualScales scales {
        50000.0, 100.0, 500000.0
    };

    const NonlinearUnknowns guess {
        340.0,
        std::log(700000.0 / scales.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto point = evaluate_operating_point(
        model, geometry, {operating, guess}, scales, options);

    if (!point.converged ||
        !point.independently_validated ||
        !point.basic_physics_valid) {
        throw std::runtime_error(
            "Operating-point prerequisite failed");
    }

    // Preserve the original solver result for comparison.
    const auto original_temperature = point.outlet_temperature;
    const auto original_pressure = point.outlet_pressure;
    const auto original_iterations = point.iterations;
    const auto original_evaluations = point.evaluations;

    const MapReferenceConditions reference {
        288.15,
        101325.0
    };

    const auto map = evaluate_compressor_map_point(
        model, geometry, point, scales, reference);

    const NonlinearUnknowns solved_unknowns {
        point.outlet_temperature,
        std::log(point.outlet_pressure / scales.pressure)
    };

    const auto full = evaluate_nonlinear_residual(
        model, geometry, operating, solved_unknowns, scales);

    const auto& inlet =
        full.full_system.system.inlet_stagnation.stagnation_state;

    const auto& outlet =
        full.full_system.system.outlet_stagnation.stagnation_state;

    check_close(
        "Inlet total temperature",
        map.inlet_total_temperature,
        inlet.T
    );

    check_close(
        "Inlet total pressure",
        map.inlet_total_pressure,
        inlet.p
    );

    check_close(
        "Outlet total temperature",
        map.outlet_total_temperature,
        outlet.T
    );

    check_close(
        "Outlet total pressure",
        map.outlet_total_pressure,
        outlet.p
    );

    check_close(
        "Total-to-total pressure ratio",
        map.total_pressure_ratio,
        outlet.p / inlet.p
    );

    check_close(
        "Static-to-static pressure ratio",
        map.static_pressure_ratio,
        point.pressure_ratio
    );

    const double tau = inlet.T / reference.temperature;
    const double delta = inlet.p / reference.pressure;

    check_close(
        "Corrected mass flow",
        map.corrected_mass_flow,
        operating.mass_flow * std::sqrt(tau) / delta
    );

    check_close(
        "Corrected speed",
        map.corrected_speed_rpm,
        operating.rotational_speed_rpm / std::sqrt(tau)
    );

    check_close(
        "Specific work",
        map.specific_work,
        point.specific_work
    );

    check_close(
        "Shaft power",
        map.shaft_power,
        point.shaft_power
    );

    if (map.operating_limits !=
        OperatingLimitStatus::NotEvaluated) {
        throw std::runtime_error(
            "Operating limits were incorrectly classified");
    }

    // The postprocessor must not modify the solver result.
    if (point.outlet_temperature != original_temperature ||
        point.outlet_pressure != original_pressure ||
        point.iterations != original_iterations ||
        point.evaluations != original_evaluations) {
        throw std::runtime_error(
            "Map postprocessing modified operating-point diagnostics");
    }

    // An unvalidated point must be rejected.
    auto invalid = point;
    invalid.independently_validated = false;

    bool rejected = false;

    try {
        (void)evaluate_compressor_map_point(
            model, geometry, invalid, scales, reference);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    if (!rejected) {
        throw std::runtime_error(
            "Unvalidated operating point was accepted");
    }

    std::cout << "PASS: " << name
              << " compressor-map integration\n";
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal;
        test_model(ideal, "IdealGas");

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        const thermomodel::PengRobinson pr(caloric);
        test_model(pr, "PengRobinson");

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
