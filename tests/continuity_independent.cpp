#include <turbocompgpu/continuity.hpp>
#include <turbocompgpu/conservation.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-12)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    const double error =
        std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) ||
        error > tolerance) {
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
void expect_failure(
    const char* name,
    Function&& function)
{
    bool failed = false;

    try {
        function();
    }
    catch (const std::exception&) {
        failed = true;
    }

    if (!failed) {
        throw std::runtime_error(name);
    }

    std::cout << "PASS: " << name << '\n';
}

InletAnnulus inlet_geometry()
{
    return {0.10, 0.20, 0.05};
}

ImpellerExitPassage exit_geometry()
{
    return {0.40, 0.025, 0.10};
}

void test_independent_mass_flows()
{
    const IndependentFlowOperatingPoint operating_point{
        5.0,
        5.5,
        20.0,
        40.0
    };

    const auto result =
        evaluate_continuity_independent(
            inlet_geometry(),
            exit_geometry(),
            operating_point);

    check_close(
        "Inlet mass flow",
        result.inlet_mass_flow,
        5.0);

    check_close(
        "Outlet mass flow",
        result.exit_mass_flow,
        5.5);

    check_close(
        "Inlet meridional velocity",
        result.inlet_meridional_velocity,
        5.0 /
            (20.0 * result.inlet_effective_area));

    check_close(
        "Outlet meridional velocity",
        result.exit_meridional_velocity,
        5.5 /
            (40.0 * result.exit_effective_area));

    const ConservationGeometry geometry{
        0.10,
        0.20
    };

    const ConservationState state{
        result.inlet_mass_flow,
        result.exit_mass_flow,
        20.0,
        200.0,
        300000.0,
        338000.0,
        190.0,
        1000.0
    };

    const auto residual =
        evaluate_conservation(geometry, state);

    check_close(
        "Independent-flow mass residual",
        residual.mass,
        0.5);

    std::cout
        << "PASS: Independent inlet/outlet mass flows\n";
}

void test_equal_flow_compatibility()
{
    const ContinuityOperatingPoint legacy{
        5.0,
        20.0,
        40.0
    };

    const IndependentFlowOperatingPoint independent{
        5.0,
        5.0,
        20.0,
        40.0
    };

    const auto original =
        evaluate_continuity(
            inlet_geometry(),
            exit_geometry(),
            legacy);

    const auto extended =
        evaluate_continuity_independent(
            inlet_geometry(),
            exit_geometry(),
            independent);

    check_close(
        "Inlet velocity compatibility",
        extended.inlet_meridional_velocity,
        original.inlet_meridional_velocity);

    check_close(
        "Outlet velocity compatibility",
        extended.exit_meridional_velocity,
        original.exit_meridional_velocity);

    check_close(
        "Inlet mass-flow compatibility",
        extended.inlet_mass_flow,
        original.inlet_mass_flow);

    check_close(
        "Outlet mass-flow compatibility",
        extended.exit_mass_flow,
        original.exit_mass_flow);

    std::cout
        << "PASS: Equal-flow compatibility\n";
}

void test_mass_flow_perturbation()
{
    IndependentFlowOperatingPoint operating_point{
        5.0,
        5.0,
        20.0,
        40.0
    };

    const auto baseline =
        evaluate_continuity_independent(
            inlet_geometry(),
            exit_geometry(),
            operating_point);

    operating_point.outlet_mass_flow += 0.5;

    const auto perturbed =
        evaluate_continuity_independent(
            inlet_geometry(),
            exit_geometry(),
            operating_point);

    check_close(
        "Inlet velocity unchanged",
        perturbed.inlet_meridional_velocity,
        baseline.inlet_meridional_velocity);

    check_close(
        "Outlet velocity perturbation",
        perturbed.exit_meridional_velocity -
            baseline.exit_meridional_velocity,
        0.5 /
            (40.0 * baseline.exit_effective_area));

    std::cout
        << "PASS: Independent outlet-flow perturbation\n";
}

void test_invalid_inputs()
{
    auto operating_point =
        IndependentFlowOperatingPoint{
            5.0, 5.0, 20.0, 40.0
        };

    operating_point.inlet_mass_flow = -1.0;

    expect_failure(
        "Negative inlet mass flow",
        [&] {
            evaluate_continuity_independent(
                inlet_geometry(),
                exit_geometry(),
                operating_point);
        });

    operating_point = {
        5.0, 5.0, 20.0, 40.0
    };

    operating_point.outlet_mass_flow =
        std::numeric_limits<double>::quiet_NaN();

    expect_failure(
        "NaN outlet mass flow",
        [&] {
            evaluate_continuity_independent(
                inlet_geometry(),
                exit_geometry(),
                operating_point);
        });

    operating_point = {
        5.0, 5.0, 20.0, 40.0
    };

    operating_point.exit_density = 0.0;

    expect_failure(
        "Zero outlet density",
        [&] {
            evaluate_continuity_independent(
                inlet_geometry(),
                exit_geometry(),
                operating_point);
        });
}

} // namespace

int main()
{
    try {
        test_independent_mass_flows();
        test_equal_flow_compatibility();
        test_mass_flow_perturbation();
        test_invalid_inputs();

        std::cout
            << "\nPASS: All independent-flow continuity tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
