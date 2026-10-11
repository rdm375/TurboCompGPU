#include <turbocompgpu/continuity.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

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

void test_geometry_and_continuity()
{
    using namespace turbocompgpu;

    const InletAnnulus inlet{
        0.10, 0.20, 0.05
    };

    const ImpellerExitPassage exit{
        0.40, 0.025, 0.10
    };

    const ContinuityOperatingPoint operating_point{
        5.0, 20.0, 40.0
    };

    const auto result =
        evaluate_continuity(
            inlet, exit, operating_point);

    const double A1g =
        std::numbers::pi / 4.0 *
        (0.20 * 0.20 - 0.10 * 0.10);

    const double A2g =
        std::numbers::pi * 0.40 * 0.025;

    const double A1 = 0.95 * A1g;
    const double A2 = 0.90 * A2g;

    const double Vm1 = 5.0 / (20.0 * A1);
    const double Vm2 = 5.0 / (40.0 * A2);

    check_close(
        "Inlet geometric area",
        result.inlet_geometric_area, A1g);

    check_close(
        "Exit geometric area",
        result.exit_geometric_area, A2g);

    check_close(
        "Inlet effective area",
        result.inlet_effective_area, A1);

    check_close(
        "Exit effective area",
        result.exit_effective_area, A2);

    check_close(
        "Inlet meridional velocity",
        result.inlet_meridional_velocity, Vm1);

    check_close(
        "Exit meridional velocity",
        result.exit_meridional_velocity, Vm2);

    check_close(
        "Inlet mass conservation",
        result.inlet_mass_flow, 5.0);

    check_close(
        "Exit mass conservation",
        result.exit_mass_flow, 5.0);

    std::cout
        << "A1 geometric = "
        << result.inlet_geometric_area << " m^2\n"
        << "A1 effective = "
        << result.inlet_effective_area << " m^2\n"
        << "A2 geometric = "
        << result.exit_geometric_area << " m^2\n"
        << "A2 effective = "
        << result.exit_effective_area << " m^2\n"
        << "Vm1 = "
        << result.inlet_meridional_velocity << " m/s\n"
        << "Vm2 = "
        << result.exit_meridional_velocity << " m/s\n";

    std::cout
        << "PASS: Geometry and continuity\n";
}

void test_impeller_adapter()
{
    using namespace turbocompgpu;

    const auto continuity =
        evaluate_continuity(
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {5.0, 20.0, 40.0});

    const ImpellerGeometry geometry{
        0.40, 45.0, 0.15
    };

    const ImpellerOperatingPoint operating_point{
        12000.0, 0.0, 0.9, 0.0, 20.0
    };

    const auto derived =
        evaluate_impeller_from_continuity(
            geometry, operating_point, continuity);

    const auto direct =
        evaluate_impeller(
            geometry,
            {
                12000.0,
                continuity.exit_meridional_velocity,
                0.9,
                continuity.inlet_meridional_velocity,
                20.0
            });

    check_close(
        "Adapter Euler work",
        derived.euler_specific_work,
        direct.euler_specific_work);

    check_close(
        "Adapter inlet relative speed",
        derived.inlet_relative_speed,
        direct.inlet_relative_speed);

    check_close(
        "Adapter exit relative speed",
        derived.relative_exit_speed,
        direct.relative_exit_speed);

    std::cout
        << "Euler work = "
        << derived.euler_specific_work
        << " J/kg\n";

    std::cout
        << "PASS: Impeller continuity adapter\n";
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

void test_invalid_inputs()
{
    using namespace turbocompgpu;

    expect_failure("Invalid inlet annulus", [] {
        inlet_annulus_area({0.20, 0.10, 0.0});
    });

    expect_failure("Zero exit width", [] {
        impeller_exit_area({0.40, 0.0, 0.0});
    });

    expect_failure("Negative blockage", [] {
        effective_flow_area(1.0, -0.1);
    });

    expect_failure("Complete blockage", [] {
        effective_flow_area(1.0, 1.0);
    });

    expect_failure("Negative density", [] {
        meridional_velocity(5.0, -20.0, 0.1);
    });

    expect_failure("Negative mass flow", [] {
        meridional_velocity(-5.0, 20.0, 0.1);
    });
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(17);

        test_geometry_and_continuity();
        test_impeller_adapter();
        test_invalid_inputs();

        std::cout
            << "PASS: All continuity tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
