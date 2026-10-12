#include <turbocompgpu/rotating_throat.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void close(
    double actual,
    double expected,
    double relative_tolerance,
    const std::string& message)
{
    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <=
            relative_tolerance *
            std::max(1.0, std::abs(expected)),
        message);
}

template<class Model>
void test_stationary_limit(
    const Model& model,
    const std::string& name)
{
    constexpr double T0 = 300.0;
    constexpr double p0 = 500000.0;
    constexpr double area = 0.002;

    const auto stationary =
        turbocompgpu::evaluate_choking_mass_flux(
            model, model, T0, p0);

    const turbocompgpu::RotatingThroatInputs input{
        T0,
        p0,
        0.15,
        0.0,
        0.0,
        area
    };

    const auto rotating =
        turbocompgpu::evaluate_rotating_throat(
            model, model, input);

    close(
        rotating.critical.pressure,
        stationary.pressure,
        1e-12,
        name + ": stationary pressure limit");

    close(
        rotating.critical.critical_mass_flux,
        stationary.critical_mass_flux,
        1e-12,
        name + ": stationary flux limit");

    close(
        rotating.critical_mass_flow,
        area * stationary.critical_mass_flux,
        1e-12,
        name + ": stationary mass flow");

    close(
        rotating.rothalpy,
        model.evaluate_Tp(T0, p0).value.h,
        1e-12,
        name + ": stationary rothalpy");

    std::cout
        << "PASS: " << name
        << " stationary-throat limit\n";
}

template<class Model>
void test_rotating_energy(
    const Model& model,
    const std::string& name)
{
    const turbocompgpu::RotatingThroatInputs input{
        300.0,
        500000.0,
        0.15,
        12000.0,
        -35.0,
        0.002
    };

    const auto result =
        turbocompgpu::evaluate_rotating_throat(
            model, model, input);

    const double U_expected =
        12000.0 *
        (2.0 * std::acos(-1.0) / 60.0) *
        0.15;

    close(
        result.blade_speed,
        U_expected,
        1e-12,
        name + ": blade speed");

    close(
        result.normal_relative_velocity,
        result.critical.sound_speed,
        1e-12,
        name + ": normal sonic condition");

    const double h_scale = std::max({
        1.0,
        std::abs(result.rothalpy),
        std::abs(result.critical.enthalpy),
        result.critical.sound_speed *
            result.critical.sound_speed
    });

    require(
        std::abs(result.rotating_energy_residual) <
            1e-9 * h_scale,
        name + ": rothalpy conservation");

    close(
        result.critical_mass_flow,
        input.effective_area *
            result.critical.critical_mass_flux,
        1e-12,
        name + ": mass flow");

    std::cout
        << std::setprecision(12)
        << name << ": U=" << result.blade_speed
        << " m/s, Wn*="
        << result.normal_relative_velocity
        << " m/s, W*="
        << result.total_relative_velocity
        << " m/s, mdot*="
        << result.critical_mass_flow
        << " kg/s, rothalpy residual="
        << result.rotating_energy_residual
        << " J/kg\n";
}

void test_invalid_geometry()
{
    thermomodel::IdealGas model;

    turbocompgpu::RotatingThroatInputs input{
        300.0, 500000.0,
        0.15, 12000.0, 0.0, 0.0
    };

    bool rejected = false;

    try {
        (void)turbocompgpu::evaluate_rotating_throat(
            model, model, input);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Zero throat area was not rejected");

    std::cout << "PASS: invalid throat geometry rejected\n";
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas air;

        thermomodel::PRParameters pr_parameters;

        thermomodel::IdealGasParameters methane_parameters;
        methane_parameters.R =
            8.31446261815324 / pr_parameters.molar_mass;
        methane_parameters.cp = 2200.0;
        methane_parameters.Tref = pr_parameters.Tref;
        methane_parameters.pref = pr_parameters.pref;

        thermomodel::IdealGas methane_ideal(
            methane_parameters);

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        pr_parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson methane_pr(
            caloric, pr_parameters);

        test_stationary_limit(air, "AirLikeIdeal");
        test_stationary_limit(methane_ideal, "MethaneIdeal");
        test_stationary_limit(methane_pr, "MethanePR");

        test_rotating_energy(air, "AirLikeIdeal");
        test_rotating_energy(methane_ideal, "MethaneIdeal");
        test_rotating_energy(methane_pr, "MethanePR");

        test_invalid_geometry();

        std::cout
            << "PASS: All rotating-throat tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
