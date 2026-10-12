#include <turbocompgpu/throat_incidence.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(double a, double b, double tolerance, const char* message)
{
    require(
        std::isfinite(a) &&
        std::abs(a - b) <=
            tolerance * std::max(1.0, std::abs(b)),
        message);
}

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

void test_independent_velocity_triangle()
{
    const ThroatIncidenceInputs input{
        5.0,
        2.0,
        9000.0,
        20.0,
        45.0
    };

    const auto result =
        evaluate_throat_incidence(geometry, input);

    const double area =
        std::numbers::pi * 0.25 *
        (0.20 * 0.20 - 0.10 * 0.10) * 0.95;

    const double Vm = 2.0 / (5.0 * area);

    const double U =
        (2.0 * std::numbers::pi * 9000.0 / 60.0) *
        (0.15 * 0.5);

    const double Wtheta = 20.0 - U;

    const double beta = std::numbers::pi / 4.0;

    const double Wn =
        Vm * std::sin(beta) -
        Wtheta * std::cos(beta);

    const double Wtau =
        Vm * std::cos(beta) +
        Wtheta * std::sin(beta);

    close(result.inlet_area, area, 1e-12, "Inlet area");
    close(result.inlet_meridional_velocity, Vm, 1e-12, "Inlet Vm");
    close(result.inlet_blade_speed, U, 1e-12, "Inlet U");
    close(result.inlet_relative_velocity.circumferential,
          Wtheta, 1e-12, "Inlet relative whirl");
    close(result.transported_local_velocity.normal,
          Wn, 1e-12, "Projected normal velocity");
    close(result.transported_local_velocity.tangential,
          Wtau, 1e-12, "Projected tangential velocity");

    close(
        relative_speed_squared(
            result.transported_relative_velocity),
        relative_speed_squared(
            result.transported_local_velocity),
        1e-12,
        "Kinetic energy invariant");

    std::cout << "PASS: independent velocity triangle\n";
}

void test_zero_incidence()
{
    const double rpm = 9000.0;
    const double density = 5.0;
    const double mdot = 2.0;

    const double area =
        std::numbers::pi * 0.25 *
        (0.20 * 0.20 - 0.10 * 0.10) * 0.95;

    const double Vm = mdot / (density * area);
    const double U =
        (2.0 * std::numbers::pi * rpm / 60.0) * 0.075;

    // For beta=45 degrees, Wtau=(Vm+Wtheta)/sqrt(2).
    // Hence Wtheta=-Vm gives zero incidence.
    const double prewhirl = U - Vm;

    const auto result =
        evaluate_throat_incidence(
            geometry,
            {density, mdot, rpm, prewhirl, 45.0});

    close(result.transported_local_velocity.tangential,
          0.0, 1e-12, "Zero tangential incidence");
    close(result.incidence_angle_deg,
          0.0, 1e-12, "Zero incidence angle");

    std::cout << "PASS: analytical zero-incidence locus\n";
}

void test_mass_flow_dependence()
{
    const ThroatIncidenceInputs input{
        5.0, 2.0, 9000.0, 20.0, 45.0
    };

    auto larger = input;
    larger.mass_flow *= 2.0;

    const auto a = evaluate_throat_incidence(geometry, input);
    const auto b = evaluate_throat_incidence(geometry, larger);

    close(b.inlet_meridional_velocity,
          2.0 * a.inlet_meridional_velocity,
          1e-12, "Mass-flow scaling");

    require(
        std::abs(
            b.transported_local_velocity.tangential -
            a.transported_local_velocity.tangential) > 1.0,
        "Incidence must depend on mass flow");

    std::cout << "PASS: mass-flow-dependent incidence\n";
}

void test_invalid_inputs()
{
    bool rejected = false;

    try {
        (void)evaluate_throat_incidence(
            geometry, {-1.0, 2.0, 9000.0, 20.0, 45.0});
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Negative density must be rejected");

    std::cout << "PASS: invalid input rejection\n";
}

} // namespace

int main()
{
    try {
        test_independent_velocity_triangle();
        test_zero_incidence();
        test_mass_flow_dependence();
        test_invalid_inputs();

        std::cout << "PASS: throat-incidence tests\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
