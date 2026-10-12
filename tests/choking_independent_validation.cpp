#include <turbocompgpu/choking_mass_flux.hpp>

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
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void check_relative(
    double actual,
    double expected,
    double tolerance,
    const std::string& message)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <= tolerance * scale,
        message);
}

// Independent fixed-pressure entropy inversion.
//
// This deliberately uses only Model::evaluate_Tp().
// It does not call derivatives_Tp() or the production
// sound-speed implementation.

thermomodel::State invert_entropy(
    const thermomodel::Model& model,
    double pressure,
    double target_entropy)
{
    const auto limits = model.metadata();

    double lo = limits.Tmin;
    double hi = limits.Tmax;

    const auto low = model.evaluate_Tp(lo, pressure);
    const auto high = model.evaluate_Tp(hi, pressure);

    require(
        static_cast<bool>(low) &&
        static_cast<bool>(high),
        "Independent entropy inversion endpoints invalid");

    const double f_lo = low.value.s - target_entropy;
    const double f_hi = high.value.s - target_entropy;

    require(
        f_lo <= 0.0 && f_hi >= 0.0,
        "Independent entropy inversion not bracketed");

    for (int iteration = 0; iteration < 160; ++iteration) {
        const double mid = 0.5 * (lo + hi);

        const auto evaluation =
            model.evaluate_Tp(mid, pressure);

        require(
            static_cast<bool>(evaluation),
            "Independent entropy evaluation failed");

        const double f =
            evaluation.value.s - target_entropy;

        if (std::abs(f) <= 1.0e-11) {
            return evaluation.value;
        }

        if (f < 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }

        if (hi - lo <=
            1.0e-13 * std::max(1.0, mid)) {
            const auto final_state =
                model.evaluate_Tp(0.5 * (lo + hi), pressure);

            require(
                static_cast<bool>(final_state),
                "Independent entropy inversion final evaluation failed");

            return final_state.value;
        }
    }

    throw std::runtime_error(
        "Independent entropy inversion failed to converge");
}

double independent_sound_speed_squared(
    const thermomodel::Model& model,
    double temperature,
    double pressure,
    double relative_step)
{
    const auto center =
        model.evaluate_Tp(temperature, pressure);

    require(
        static_cast<bool>(center),
        "Independent sound-speed reference state invalid");

    const double p_plus =
        pressure * (1.0 + relative_step);

    const double p_minus =
        pressure * (1.0 - relative_step);

    const auto plus =
        invert_entropy(model, p_plus, center.value.s);

    const auto minus =
        invert_entropy(model, p_minus, center.value.s);

    const double density_difference =
        plus.rho - minus.rho;

    require(
        density_difference > 0.0,
        "Nonpositive isentropic density increment");

    return (p_plus - p_minus) / density_difference;
}

void test_peng_robinson_sound_speed()
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

    constexpr std::array<double, 3> steps{
        1.0e-3, 1.0e-4, 1.0e-5
    };

    std::cout << std::setprecision(12);

    for (const auto& point : points) {
        const auto analytical =
            turbocompgpu::evaluate_isentropic_sound_speed(
                model, point.T, point.p);

        double best_error = 1.0;

        for (double step : steps) {
            const double finite_difference =
                independent_sound_speed_squared(
                    model, point.T, point.p, step);

            const double error =
                std::abs(
                    finite_difference -
                    analytical.sound_speed_squared
                ) / analytical.sound_speed_squared;

            best_error = std::min(best_error, error);
        }

        require(
            best_error < 1.0e-5,
            "Peng-Robinson independent sound-speed mismatch");

        std::cout
            << "PR sound-speed oracle: T=" << point.T
            << " K, p=" << point.p
            << " Pa, a=" << analytical.sound_speed
            << " m/s, best relative a2 error="
            << best_error << '\n';
    }

    std::cout
        << "PASS: independent Peng-Robinson "
        << "sound-speed validation\n";
}

double independent_mass_flux(
    const thermomodel::Model& model,
    double stagnation_enthalpy,
    double stagnation_entropy,
    double pressure)
{
    const auto state =
        invert_entropy(
            model, pressure, stagnation_entropy);

    const double kinetic_energy =
        stagnation_enthalpy - state.h;

    require(
        kinetic_energy > 0.0,
        "Nonpositive kinetic energy in mass-flux oracle");

    return state.rho *
           std::sqrt(2.0 * kinetic_energy);
}

template<class Model>
void test_mass_flux_maximum(
    const Model& model,
    double stagnation_temperature,
    double stagnation_pressure,
    const std::string& label)
{
    const auto stagnation =
        model.evaluate_Tp(
            stagnation_temperature,
            stagnation_pressure);

    require(
        static_cast<bool>(stagnation),
        "Invalid stagnation reference state");

    const auto critical =
        turbocompgpu::evaluate_choking_mass_flux(
            model,
            model,
            stagnation_temperature,
            stagnation_pressure);

    const double h0 = stagnation.value.h;
    const double s0 = stagnation.value.s;

    const double p_star = critical.pressure;

    const double G_center =
        independent_mass_flux(
            model, h0, s0, p_star);

    check_relative(
        G_center,
        critical.critical_mass_flux,
        1.0e-7,
        label + ": critical mass-flux mismatch");

    // Use a relatively wide perturbation to make the
    // difference in flux larger than floating-point noise.
    constexpr double pressure_fraction = 0.01;

    const double G_lower =
        independent_mass_flux(
            model,
            h0,
            s0,
            p_star * (1.0 - pressure_fraction));

    const double G_upper =
        independent_mass_flux(
            model,
            h0,
            s0,
            p_star * (1.0 + pressure_fraction));

    require(
        G_center > G_lower,
        label + ": lower-pressure flux exceeds critical flux");

    require(
        G_center > G_upper,
        label + ": upper-pressure flux exceeds critical flux");

    check_relative(
        critical.energy_residual,
        0.0,
        1.0e-4,
        label + ": sonic energy residual");

    check_relative(
        critical.entropy_residual,
        0.0,
        1.0e-8,
        label + ": sonic entropy residual");

    std::cout
        << label
        << ": G(-1%)=" << G_lower
        << ", G*=" << G_center
        << ", G(+1%)=" << G_upper
        << " kg/(m2 s)\n";
}

void test_choking_maxima()
{
    thermomodel::IdealGas ideal;

    test_mass_flux_maximum(
        ideal,
        300.0,
        500000.0,
        "AirLikeIdeal");

    thermomodel::IdealGasParameters methane_parameters;

    constexpr double universal_R =
        8.31446261815324;

    thermomodel::PRParameters pr_parameters;

    methane_parameters.R =
        universal_R / pr_parameters.molar_mass;

    methane_parameters.cp = 2200.0;
    methane_parameters.Tref = pr_parameters.Tref;
    methane_parameters.pref = pr_parameters.pref;

    thermomodel::IdealGas methane_ideal(
        methane_parameters);

    test_mass_flux_maximum(
        methane_ideal,
        300.0,
        500000.0,
        "MethaneIdeal");

    auto caloric =
        std::make_shared<thermomodel::ConstantCp>(2200.0);

    pr_parameters.root_policy =
        thermomodel::PRRootPolicy::Vapor;

    thermomodel::PengRobinson methane_pr(
        caloric, pr_parameters);

    test_mass_flux_maximum(
        methane_pr,
        300.0,
        500000.0,
        "MethanePR");

    std::cout
        << "PASS: independent choking mass-flux maxima\n";
}

} // namespace

int main()
{
    try {
        test_peng_robinson_sound_speed();
        test_choking_maxima();

        std::cout
            << "PASS: All independent choking validations\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
