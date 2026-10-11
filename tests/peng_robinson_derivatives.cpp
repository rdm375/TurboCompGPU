#include <turbocompgpu/thermodynamic_solver.hpp>

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

struct TestPoint {
    double T;
    double pressure;
};

double relative_error(double actual, double expected)
{
    return std::abs(actual - expected) /
           std::max(std::abs(expected), 1.0);
}

void validate_point(
    const thermomodel::PengRobinson& model,
    const TestPoint& point)
{
    const auto center = turbocompgpu::evaluate(
        model, point.T, point.pressure);

    if (!std::isfinite(center.cp) || center.cp <= 0.0) {
        throw std::runtime_error("Invalid analytical cp");
    }

    constexpr std::array<double, 5> steps{
        1.0e-1,
        1.0e-2,
        1.0e-3,
        1.0e-4,
        1.0e-5
    };

    std::cout
        << "\nT = " << point.T
        << " K, P = " << point.pressure / 1.0e6
        << " MPa\n"
        << "cp = " << center.cp << " J/(kg K)\n";

    std::cout
        << std::setw(12) << "dT"
        << std::setw(20) << "cp from dh/dT"
        << std::setw(20) << "cp from T ds/dT"
        << std::setw(20) << "h relative error"
        << std::setw(20) << "s relative error"
        << '\n';

    double best_h_error = 1.0;
    double best_s_error = 1.0;

    for (const double step : steps) {
        const auto plus = turbocompgpu::evaluate(
            model, point.T + step, point.pressure);

        const auto minus = turbocompgpu::evaluate(
            model, point.T - step, point.pressure);

        const double cp_h =
            (plus.h - minus.h) / (2.0 * step);

        const double cp_s =
            point.T * (plus.s - minus.s) / (2.0 * step);

        const double error_h =
            relative_error(cp_h, center.cp);

        const double error_s =
            relative_error(cp_s, center.cp);

        best_h_error = std::min(best_h_error, error_h);
        best_s_error = std::min(best_s_error, error_s);

        std::cout
            << std::setw(12) << step
            << std::setw(20) << cp_h
            << std::setw(20) << cp_s
            << std::setw(20) << error_h
            << std::setw(20) << error_s
            << '\n';
    }

    constexpr double tolerance = 1.0e-5;

    if (best_h_error > tolerance ||
        best_s_error > tolerance) {
        throw std::runtime_error(
            "Peng-Robinson heat-capacity derivative mismatch");
    }

    std::cout
        << "PASS: derivative consistency\n";
}

} // namespace

int main()
{
    try {
        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson model(
            caloric, parameters);

        constexpr std::array<TestPoint, 6> points{{
            {250.0, 1.0e6},
            {300.0, 5.0e6},
            {350.0, 5.0e6},
            {350.0, 10.0e6},
            {400.0, 10.0e6},
            {500.0, 15.0e6}
        }};

        std::cout << std::setprecision(12);

        for (const auto& point : points) {
            validate_point(model, point);
        }

        std::cout
            << "\nPASS: All Peng-Robinson derivative tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
