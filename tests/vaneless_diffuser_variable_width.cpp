#include <turbocompgpu/vaneless_diffuser_variable_width.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

double relative_error(double a, double b)
{
    return std::abs(a - b) /
        std::max({1.0, std::abs(a), std::abs(b)});
}

template<class Model>
void run_model(const Model& model, const char* name)
{
    constexpr double r2 = 0.20;
    constexpr double r3 = 0.30;
    constexpr double b2 = 0.025;
    constexpr double b3 = 0.020;

    const auto state = evaluate(model, 340.0, 700000.0);

    constexpr double u2 = 60.0;
    constexpr double v2 = 90.0;

    const double mass =
        state.rho *
        (2.0 * std::numbers::pi * r2 * b2) *
        u2;

    const VanelessDiffuserInlet inlet{
        state, u2, v2, mass
    };

    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = 0.01;
    options.radial_steps = 512;

    const VanelessDiffuserGeometry constant_geometry{
        r2, r3, b2, 0.0
    };

    const VariableWidthDiffuserGeometry constant_profile{
        r2, r3,
        [](double) { return b2; },
        [](double) { return 0.0; }
    };

    const auto reference =
        evaluate_frictional_vaneless_diffuser_analytical(
            model, model, constant_geometry, inlet, options);

    const auto constant =
        evaluate_variable_width_vaneless_diffuser(
            model, model, constant_profile, inlet, options);

    require(
        relative_error(
            constant.outlet_static.p,
            reference.outlet_static.p) < 1.0e-12,
        "Constant-width pressure equivalence failed");

    require(
        relative_error(
            constant.outlet_static.T,
            reference.outlet_static.T) < 1.0e-12,
        "Constant-width temperature equivalence failed");

    require(
        relative_error(
            constant.outlet_tangential_velocity,
            reference.outlet_tangential_velocity) < 1.0e-12,
        "Constant-width whirl equivalence failed");

    const VariableWidthDiffuserGeometry tapered{
        r2, r3,
        [](double r) {
            return b2 +
                (b3 - b2) * (r - r2) / (r3 - r2);
        },
        [](double) {
            return (b3 - b2) / (r3 - r2);
        }
    };

    const auto result =
        evaluate_variable_width_vaneless_diffuser(
            model, model, tapered, inlet, options);

    const auto finite_difference =
        evaluate_variable_width_vaneless_diffuser_finite_difference(
            model, tapered, inlet, options);

    const double pressure_difference = relative_error(
        result.outlet_static.p,
        finite_difference.outlet_static.p);

    const double temperature_difference = relative_error(
        result.outlet_static.T,
        finite_difference.outlet_static.T);

    const double whirl_difference = relative_error(
        result.outlet_tangential_velocity,
        finite_difference.outlet_tangential_velocity);

    require(
        pressure_difference < 1.0e-7,
        "Variable-width analytical/FD pressure mismatch");

    require(
        temperature_difference < 1.0e-7,
        "Variable-width analytical/FD temperature mismatch");

    require(
        whirl_difference < 1.0e-7,
        "Variable-width analytical/FD whirl mismatch");

    std::cout
        << name
        << " analytical_vs_fd"
        << " pressure_error=" << pressure_difference
        << " temperature_error=" << temperature_difference
        << " whirl_error=" << whirl_difference
        << '\n';

    require(
        std::abs(result.mass_residual) <
            1.0e-10 * std::max(1.0, mass),
        "Variable-width mass conservation failed");

    require(
        std::abs(result.stagnation_enthalpy_residual) <
            1.0e-7,
        "Variable-width energy conservation failed");

    require(
        std::abs(result.entropy_balance_residual) <
            1.0e-7,
        "Variable-width entropy balance failed");

    require(
        std::abs(result.angular_momentum_balance_residual) <
            1.0e-8,
        "Variable-width angular momentum balance failed");

    require(
        result.entropy_change > 0.0,
        "Variable-width entropy production must be positive");

    std::cout
        << name
        << " constant_width=PASS"
        << " tapered_p3=" << result.outlet_static.p
        << " tapered_T3=" << result.outlet_static.T
        << " energy_residual="
        << result.stagnation_enthalpy_residual
        << " entropy_residual="
        << result.entropy_balance_residual
        << '\n';
}

} // namespace

int main()
{
    try {
        std::cout << std::scientific
                  << std::setprecision(12);

        thermomodel::IdealGas ideal;
        run_model(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(caloric, parameters);
        run_model(pr, "PengRobinson");

        std::cout
            << "Variable-width diffuser initial verification: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
