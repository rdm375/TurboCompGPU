#include <turbocompgpu/impeller.hpp>

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

void test_radial_blades()
{
    const turbocompgpu::ImpellerGeometry geometry{
        0.4, 90.0
    };

    const turbocompgpu::ImpellerOperatingPoint operating_point{
        12000.0, 100.0, 1.0
    };

    const auto result =
        turbocompgpu::evaluate_impeller(
            geometry, operating_point);

    const double U2 =
        std::numbers::pi * 0.4 * 12000.0 / 60.0;

    check_close("Radial tip speed", result.tip_speed, U2);

    check_close(
        "Radial ideal whirl",
        result.ideal_exit_whirl,
        U2);

    check_close(
        "Radial actual whirl",
        result.actual_exit_whirl,
        U2);

    check_close(
        "Radial Euler work",
        result.euler_specific_work,
        U2 * U2);

    check_close(
        "Radial relative speed",
        result.relative_exit_speed,
        100.0);

    check_close(
        "Radial relative angle",
        result.relative_exit_angle_deg,
        90.0);

    std::cout << "PASS: Radial-blade analytical case\n";
}

void test_backward_curved()
{
    const turbocompgpu::ImpellerGeometry geometry{
        0.4, 45.0
    };

    const turbocompgpu::ImpellerOperatingPoint operating_point{
        12000.0, 100.0, 0.9
    };

    const auto result =
        turbocompgpu::evaluate_impeller(
            geometry, operating_point);

    const double U2 =
        std::numbers::pi * 0.4 * 12000.0 / 60.0;

    const double expected_whirl =
        0.9 * (U2 - 100.0);

    check_close(
        "Backward-curved ideal whirl",
        result.ideal_exit_whirl,
        U2 - 100.0);

    check_close(
        "Backward-curved actual whirl",
        result.actual_exit_whirl,
        expected_whirl);

    check_close(
        "Backward-curved Euler work",
        result.euler_specific_work,
        U2 * expected_whirl);

    const double expected_relative_angle =
        std::atan2(
            100.0,
            U2 - expected_whirl) *
        180.0 / std::numbers::pi;

    check_close(
        "Backward-curved relative angle",
        result.relative_exit_angle_deg,
        expected_relative_angle);

    std::cout
        << "PASS: Backward-curved analytical case\n";
}

void test_slip_scaling()
{
    const turbocompgpu::ImpellerGeometry geometry{
        0.4, 60.0
    };

    const auto no_slip =
        turbocompgpu::evaluate_impeller(
            geometry, {12000.0, 100.0, 1.0});

    const auto slip =
        turbocompgpu::evaluate_impeller(
            geometry, {12000.0, 100.0, 0.8});

    check_close(
        "Slip factor scales whirl",
        slip.actual_exit_whirl,
        0.8 * no_slip.actual_exit_whirl);

    check_close(
        "Slip factor scales Euler work",
        slip.euler_specific_work,
        0.8 * no_slip.euler_specific_work);

    std::cout << "PASS: Slip-factor scaling\n";
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

    expect_failure("Negative diameter", [] {
        evaluate_impeller(
            {-0.4, 45.0},
            {12000.0, 100.0, 0.9});
    });

    expect_failure("Zero blade angle", [] {
        evaluate_impeller(
            {0.4, 0.0},
            {12000.0, 100.0, 0.9});
    });

    expect_failure("Slip factor above unity", [] {
        evaluate_impeller(
            {0.4, 45.0},
            {12000.0, 100.0, 1.1});
    });

    expect_failure("Negative speed", [] {
        evaluate_impeller(
            {0.4, 45.0},
            {-12000.0, 100.0, 0.9});
    });

    expect_failure("Negative meridional velocity", [] {
        evaluate_impeller(
            {0.4, 45.0},
            {12000.0, -100.0, 0.9});
    });
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(17);

        test_radial_blades();
        test_backward_curved();
        test_slip_scaling();
        test_invalid_inputs();

        std::cout
            << "PASS: All impeller aerodynamic tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
