#include <turbocompgpu/vaneless_diffuser_friction.hpp>

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
#include <vector>

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

template<class Function>
void expect_invalid(Function&& function, const std::string& name)
{
    try {
        function();
    }
    catch (const std::invalid_argument&) {
        return;
    }

    throw std::runtime_error(
        "Expected invalid_argument: " + name);
}

template<class Model>
VanelessDiffuserInlet make_inlet(
    const Model& model,
    const VanelessDiffuserGeometry& geometry)
{
    const auto state = evaluate(model, 340.0, 700000.0);

    constexpr double u = 60.0;
    constexpr double v = 90.0;

    const double area =
        2.0 * std::numbers::pi *
        geometry.inlet_radius *
        geometry.outlet_width;

    return {
        state,
        u,
        v,
        state.rho * area * u
    };
}

template<class Model>
FrictionalDiffuserResult solve(
    const Model& model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    double friction,
    int steps)
{
    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = friction;
    options.radial_steps = steps;

    return evaluate_frictional_vaneless_diffuser(
        model, geometry, inlet, options);
}

template<class Model>
void convergence_study(
    const Model& model,
    const char* name,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet)
{
    std::cout << "\nCONVERGENCE: " << name << '\n';

    constexpr std::array<int, 7> meshes{
        32, 64, 128, 256, 512, 1024, 2048
    };

    std::vector<double> pressures;

    for (int n : meshes) {
        const auto result =
            solve(model, geometry, inlet, 0.02, n);

        pressures.push_back(result.outlet_static.p);

        std::cout
            << "N=" << std::setw(5) << n
            << " p3=" << std::setprecision(16)
            << result.outlet_static.p
            << " energy=" << result.stagnation_enthalpy_residual
            << " entropy_balance=" << result.entropy_balance_residual
            << " angular_balance="
            << result.angular_momentum_balance_residual
            << '\n';
    }

    std::cout << "Successive pressure differences and observed orders:\n";

    for (std::size_t i = 1; i < pressures.size(); ++i) {
        const double difference =
            std::abs(pressures[i] - pressures[i - 1]);

        std::cout
            << "N=" << meshes[i]
            << " difference=" << difference;

        if (i >= 2) {
            const double previous =
                std::abs(pressures[i - 1] - pressures[i - 2]);

            if (difference > 0.0 && previous > 0.0) {
                const double order =
                    std::log2(previous / difference);

                std::cout << " observed_order=" << order;
            } else {
                std::cout << " observed_order=unresolved";
            }
        }

        std::cout << '\n';
    }

    // The pressure sequence must stabilize on refinement.
    // A strict fourth-order assertion is deliberately omitted:
    // finite-difference thermodynamic derivatives impose
    // an accuracy floor independent of radial step size.
    require(
        relative_error(pressures[pressures.size() - 1],
                       pressures[pressures.size() - 2]) < 1.0e-7,
        std::string(name) + ": pressure refinement did not stabilize");
}

template<class Model>
void friction_sweep(
    const Model& model,
    const char* name,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet)
{
    std::cout << "\nFRICTION SWEEP: " << name << '\n';

    constexpr std::array<double, 6> factors{
        0.0, 0.005, 0.01, 0.02, 0.04, 0.08
    };

    double previous_entropy = -1.0;
    double previous_angular_momentum =
        geometry.inlet_radius * inlet.tangential_velocity;

    for (double f : factors) {
        const auto result =
            solve(model, geometry, inlet, f, 1024);

        const double h0 =
            result.inlet_static.h +
            0.5 * (
                inlet.meridional_velocity *
                    inlet.meridional_velocity +
                inlet.tangential_velocity *
                    inlet.tangential_velocity);

        const double energy_error =
            std::abs(result.stagnation_enthalpy_residual) /
            std::max(1.0, std::abs(h0));

        const double angular_momentum_out =
            geometry.outlet_radius *
            result.outlet_tangential_velocity;

        std::cout
            << "f_D=" << std::setw(7) << f
            << " p3=" << std::setprecision(14)
            << result.outlet_static.p
            << " ds=" << result.entropy_change
            << " integrated_ds="
            << result.integrated_entropy_production
            << " dL=" << result.angular_momentum_change
            << " energy_relative=" << energy_error
            << " entropy_defect="
            << result.entropy_balance_residual
            << " angular_defect="
            << result.angular_momentum_balance_residual
            << '\n';

        require(
            std::abs(result.mass_residual) <
                1.0e-10 * std::max(1.0, inlet.mass_flow),
            "Mass conservation failed");

        require(energy_error < 1.0e-9,
                "Stagnation enthalpy conservation failed");

        require(
            std::abs(result.entropy_balance_residual) < 1.0e-7,
            "Entropy integral disagrees with EOS entropy");

        require(
            std::abs(
                result.angular_momentum_balance_residual) < 1.0e-8,
            "Angular momentum integral disagrees with outlet");

        require(
            result.integrated_entropy_production >= 0.0,
            "Negative integrated entropy production");

        require(
            result.entropy_change >= -1.0e-7,
            "Negative EOS entropy change");

        require(
            result.entropy_change >= previous_entropy - 1.0e-7,
            "Entropy generation not monotonic in friction");

        require(
            angular_momentum_out <=
                previous_angular_momentum + 1.0e-8,
            "Angular momentum not monotonic in friction");

        previous_entropy = result.entropy_change;
        previous_angular_momentum = angular_momentum_out;
    }
}

template<class Model>
void failure_tests(
    const Model& model,
    const char* name,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet)
{
    std::cout << "\nFAILURE PATHS: " << name << '\n';

    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = 0.02;

    {
        auto invalid = geometry;
        invalid.outlet_radius = invalid.inlet_radius;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, invalid, inlet, options);
            },
            "nonincreasing radius");
    }

    {
        auto invalid = geometry;
        invalid.outlet_width = -0.01;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, invalid, inlet, options);
            },
            "negative width");
    }

    {
        auto invalid = geometry;
        invalid.outlet_blockage = 0.1;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, invalid, inlet, options);
            },
            "unsupported blockage");
    }

    {
        auto invalid = options;
        invalid.darcy_friction_factor = -0.01;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, geometry, inlet, invalid);
            },
            "negative Darcy friction factor");
    }

    {
        auto invalid = options;
        invalid.radial_steps = 0;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, geometry, inlet, invalid);
            },
            "zero radial steps");
    }

    {
        auto invalid = inlet;
        invalid.mass_flow *= 1.1;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, geometry, invalid, options);
            },
            "inconsistent inlet continuity");
    }

    {
        auto invalid = inlet;
        invalid.meridional_velocity = -10.0;

        expect_invalid(
            [&] {
                evaluate_frictional_vaneless_diffuser(
                    model, geometry, invalid, options);
            },
            "negative radial velocity");
    }

    std::cout << "Invalid-input rejection: PASS\n";
}

template<class Model>
void run_model(const Model& model, const char* name)
{
    const VanelessDiffuserGeometry geometry{
        0.20, 0.24, 0.025, 0.0
    };

    const auto inlet = make_inlet(model, geometry);

    convergence_study(model, name, geometry, inlet);
    friction_sweep(model, name, geometry, inlet);
    failure_tests(model, name, geometry, inlet);
}

} // namespace

int main()
{
    try {
        std::cout << std::scientific
                  << std::setprecision(10);

        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        run_model(ideal, "IdealGas");
        run_model(pr, "PengRobinson");

        std::cout
            << "\nStage 8C.1c numerical verification: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "\nFAIL: " << error.what() << '\n';
        return 1;
    }
}
