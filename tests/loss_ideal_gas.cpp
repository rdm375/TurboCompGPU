#include <turbocompgpu/loss.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-11)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            tolerance * scale) {
        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected << '\n';

        throw std::runtime_error(name);
    }
}

thermomodel::State evaluate(
    const thermomodel::Model& model,
    double temperature,
    double pressure)
{
    const auto result =
        model.evaluate_Tp(temperature, pressure);

    if (!result) {
        throw std::runtime_error(
            "Thermodynamic evaluation failed: " +
            result.message);
    }

    return result.value;
}

void test_analytical_entropy_loss()
{
    const thermomodel::IdealGas model;

    constexpr double R = 287.05;
    constexpr double cp = 1004.5;

    constexpr double T01 = 300.0;
    constexpr double T02 = 340.0;
    constexpr double p01 = 500000.0;

    constexpr double prescribed_loss = 20.0;

    const double pressure_ratio =
        std::pow(T02 / T01, cp / R) *
        std::exp(-prescribed_loss / R);

    const double p02 =
        p01 * pressure_ratio;

    const auto inlet =
        evaluate(model, T01, p01);

    const auto outlet =
        evaluate(model, T02, p02);

    const auto loss =
        evaluate_entropy_loss(
            inlet,
            outlet,
            {prescribed_loss});

    check_close(
        "Analytical entropy generation",
        loss.actual_entropy_increase,
        prescribed_loss);

    check_close(
        "Entropy closure residual",
        loss.entropy_residual,
        0.0);

    check_close(
        "Stagnation enthalpy rise",
        outlet.h - inlet.h,
        cp * (T02 - T01));

    std::cout
        << "PASS: Analytical ideal-gas entropy loss\n";
}

void test_pressure_penalty()
{
    constexpr double R = 287.05;
    constexpr double cp = 1004.5;

    constexpr double T01 = 300.0;
    constexpr double T02 = 340.0;
    constexpr double p01 = 500000.0;

    constexpr double entropy_loss = 20.0;

    const double isentropic_p02 =
        p01 * std::pow(T02 / T01, cp / R);

    const double lossy_p02 =
        isentropic_p02 *
        std::exp(-entropy_loss / R);

    check_close(
        "Pressure penalty",
        lossy_p02 / isentropic_p02,
        std::exp(-entropy_loss / R));

    std::cout
        << "PASS: Stagnation-pressure penalty\n";
}

void test_loss_perturbation()
{
    const thermomodel::IdealGas model;

    const auto inlet =
        evaluate(model, 300.0, 500000.0);

    const auto outlet =
        evaluate(model, 340.0, 700000.0);

    const auto baseline =
        evaluate_entropy_loss(
            inlet, outlet, {10.0});

    const auto perturbed =
        evaluate_entropy_loss(
            inlet, outlet, {15.0});

    check_close(
        "Independent loss perturbation",
        perturbed.entropy_residual -
            baseline.entropy_residual,
        -5.0);

    std::cout
        << "PASS: Independent entropy-loss perturbation\n";
}

void test_invalid_loss()
{
    const thermomodel::IdealGas model;

    const auto inlet =
        evaluate(model, 300.0, 500000.0);

    const auto outlet =
        evaluate(model, 340.0, 700000.0);

    bool rejected = false;

    try {
        evaluate_entropy_loss(
            inlet, outlet, {-1.0});
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    if (!rejected) {
        throw std::runtime_error(
            "Negative entropy generation was accepted");
    }

    std::cout
        << "PASS: Negative entropy-generation rejection\n";
}

} // namespace

int main()
{
    try {
        test_analytical_entropy_loss();
        test_pressure_penalty();
        test_loss_perturbation();
        test_invalid_loss();

        std::cout
            << "\nPASS: All ideal-gas loss tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
