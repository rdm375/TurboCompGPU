#include <turbocompgpu/isentropic_sound_speed.hpp>

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

void require_close(
    double actual,
    double expected,
    double tolerance,
    const char* message)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        throw std::runtime_error(message);
    }
}

void test_ideal_gas()
{
    thermomodel::IdealGas model;

    const auto& parameters = model.parameters();

    const double R = parameters.R;
    const double cp = parameters.cp;
    const double cv = cp - R;
    const double gamma = cp / cv;

    if (!(cv > 0.0 && gamma > 1.0)) {
        throw std::runtime_error(
            "Invalid ideal-gas heat capacities");
    }

    constexpr std::array<double, 4> temperatures{
        200.0, 300.0, 500.0, 1000.0
    };

    constexpr std::array<double, 3> pressures{
        1.0e5, 5.0e5, 5.0e6
    };

    for (double T : temperatures) {
        for (double p : pressures) {
            const auto result =
                turbocompgpu::evaluate_isentropic_sound_speed(
                    model, T, p);

            const double expected_a2 = gamma * R * T;
            const double expected_a = std::sqrt(expected_a2);

            require_close(
                result.sound_speed_squared,
                expected_a2,
                1.0e-11,
                "Ideal-gas sound-speed-squared mismatch");

            require_close(
                result.sound_speed,
                expected_a,
                1.0e-11,
                "Ideal-gas sound-speed mismatch");

            const auto state = model.evaluate_Tp(T, p);

            if (!state) {
                throw std::runtime_error(
                    "Ideal-gas state evaluation failed");
            }

            if (std::isfinite(state.value.sound_speed)) {
                require_close(
                    result.sound_speed,
                    state.value.sound_speed,
                    1.0e-11,
                    "Ideal-gas EOS sound-speed mismatch");
            }
        }
    }

    std::cout
        << "PASS: ideal-gas analytical sound speed (12 states)\n";
}

void test_peng_robinson()
{
    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(2200.0);

    thermomodel::PRParameters parameters;
    parameters.root_policy =
        thermomodel::PRRootPolicy::Vapor;

    thermomodel::PengRobinson model(
        caloric, parameters);

    struct Point {
        double T;
        double p;
    };

    constexpr std::array<Point, 6> points{{
        {250.0, 1.0e6},
        {300.0, 5.0e6},
        {350.0, 5.0e6},
        {350.0, 10.0e6},
        {400.0, 10.0e6},
        {500.0, 15.0e6}
    }};

    std::cout << std::setprecision(12);

    for (const auto& point : points) {
        const auto result =
            turbocompgpu::evaluate_isentropic_sound_speed(
                model, point.T, point.p);

        const auto state =
            model.evaluate_Tp(point.T, point.p);

        if (!state) {
            throw std::runtime_error(
                "Peng-Robinson state evaluation failed");
        }

        if (!(result.sound_speed > 0.0 &&
              result.sound_speed_squared > 0.0)) {
            throw std::runtime_error(
                "Invalid Peng-Robinson sound speed");
        }

        std::cout
            << "PR: T=" << point.T
            << " K, p=" << point.p
            << " Pa, rho=" << state.value.rho
            << " kg/m3, a=" << result.sound_speed
            << " m/s\n";
    }

    std::cout
        << "PASS: Peng-Robinson analytical-derivative "
        << "sound speed (6 states)\n";
}

} // namespace

int main()
{
    try {
        test_ideal_gas();
        test_peng_robinson();

        std::cout
            << "PASS: All isentropic sound-speed tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
