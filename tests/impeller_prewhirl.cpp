#include <turbocompgpu/impeller.hpp>

#include <algorithm>
#include <array>
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

void test_prewhirl()
{
    using namespace turbocompgpu;

    const ImpellerGeometry geometry{
        0.4, 45.0, 0.2
    };

    const double N = 12000.0;
    const double Vm1 = 80.0;
    const double Vm2 = 100.0;
    const double sigma = 0.9;

    const double U1 =
        std::numbers::pi * 0.2 * N / 60.0;

    const double U2 =
        std::numbers::pi * 0.4 * N / 60.0;

    const auto baseline =
        evaluate_impeller(
            geometry,
            {N, Vm2, sigma, Vm1, 0.0});

    check_close(
        "Inlet blade speed",
        baseline.inlet_blade_speed,
        U1);

    check_close(
        "Exit blade speed",
        baseline.tip_speed,
        U2);

    check_close(
        "Zero-prewhirl Euler work",
        baseline.euler_specific_work,
        baseline.exit_euler_work);

    const std::array<double, 5> prewhirl_values{
        -40.0, -20.0, 0.0, 20.0, 40.0
    };

    std::cout
        << std::setw(15) << "Vtheta1"
        << std::setw(20) << "Euler work"
        << std::setw(20) << "Delta work"
        << std::setw(20) << "Expected delta"
        << '\n';

    for (double Vtheta1 : prewhirl_values) {
        const auto result =
            evaluate_impeller(
                geometry,
                {N, Vm2, sigma, Vm1, Vtheta1});

        const double delta_work =
            result.euler_specific_work -
            baseline.euler_specific_work;

        const double expected_delta =
            -U1 * Vtheta1;

        std::cout
            << std::setw(15) << Vtheta1
            << std::setw(20) << result.euler_specific_work
            << std::setw(20) << delta_work
            << std::setw(20) << expected_delta
            << '\n';

        check_close(
            "Euler prewhirl correction",
            delta_work,
            expected_delta);

        check_close(
            "Inlet Euler contribution",
            result.inlet_euler_work,
            U1 * Vtheta1);

        check_close(
            "Inlet absolute speed",
            result.inlet_absolute_speed,
            std::hypot(Vm1, Vtheta1));

        check_close(
            "Inlet relative speed",
            result.inlet_relative_speed,
            std::hypot(Vm1, Vtheta1 - U1));

        const double expected_angle =
            std::atan2(Vm1, U1 - Vtheta1) *
            180.0 / std::numbers::pi;

        check_close(
            "Inlet relative flow angle",
            result.inlet_relative_angle_deg,
            expected_angle);
    }

    std::cout
        << "PASS: Inlet prewhirl analytical validation\n";
}

void test_zero_prewhirl_compatibility()
{
    using namespace turbocompgpu;

    const auto legacy =
        evaluate_impeller(
            {0.4, 45.0},
            {12000.0, 100.0, 0.9});

    const auto explicit_zero =
        evaluate_impeller(
            {0.4, 45.0, 0.2},
            {12000.0, 100.0, 0.9, 80.0, 0.0});

    check_close(
        "Zero-prewhirl compatibility",
        explicit_zero.euler_specific_work,
        legacy.euler_specific_work);

    std::cout
        << "PASS: Zero-prewhirl compatibility\n";
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

    expect_failure("Negative inlet diameter", [] {
        evaluate_impeller(
            {0.4, 45.0, -0.1},
            {12000.0, 100.0, 0.9, 80.0, 0.0});
    });

    expect_failure("Inlet diameter exceeds exit", [] {
        evaluate_impeller(
            {0.4, 45.0, 0.5},
            {12000.0, 100.0, 0.9, 80.0, 0.0});
    });

    expect_failure("Prewhirl without inlet diameter", [] {
        evaluate_impeller(
            {0.4, 45.0},
            {12000.0, 100.0, 0.9, 80.0, 20.0});
    });

    expect_failure("Negative inlet meridional velocity", [] {
        evaluate_impeller(
            {0.4, 45.0, 0.2},
            {12000.0, 100.0, 0.9, -80.0, 0.0});
    });
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(17);

        test_prewhirl();
        test_zero_prewhirl_compatibility();
        test_invalid_inputs();

        std::cout
            << "PASS: All inlet-prewhirl tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
