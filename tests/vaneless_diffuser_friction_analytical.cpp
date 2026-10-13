#include <turbocompgpu/vaneless_diffuser_friction.hpp>
#include <turbocompgpu/vaneless_diffuser_friction_analytical.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

double relative_error(double a, double b)
{
    return std::abs(a - b) /
        std::max({1.0, std::abs(a), std::abs(b)});
}

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template<class Model>
void verify(const Model& model, const char* name)
{
    const VanelessDiffuserGeometry geometry{
        0.20, 0.24, 0.025, 0.0
    };

    const auto state = evaluate(model, 340.0, 700000.0);

    constexpr double u = 60.0;
    constexpr double v = 90.0;

    const double area =
        2.0 * std::numbers::pi *
        geometry.inlet_radius *
        geometry.outlet_width;

    const VanelessDiffuserInlet inlet{
        state, u, v, state.rho * area * u
    };

    constexpr std::array<double, 6> factors{
        0.0, 0.005, 0.01, 0.02, 0.04, 0.08
    };

    constexpr std::array<int, 4> meshes{
        32, 64, 128, 512
    };

    double maximum_pressure_error = 0.0;
    double maximum_temperature_error = 0.0;
    double maximum_velocity_error = 0.0;

    for (double friction : factors) {
        for (int steps : meshes) {
            FrictionalDiffuserOptions options;
            options.darcy_friction_factor = friction;
            options.radial_steps = steps;

            const auto reference =
                evaluate_frictional_vaneless_diffuser(
                    model, geometry, inlet, options);

            const auto analytical =
                evaluate_frictional_vaneless_diffuser_analytical(
                    model, model, geometry, inlet, options);

            const double ep = relative_error(
                analytical.outlet_static.p,
                reference.outlet_static.p);

            const double eT = relative_error(
                analytical.outlet_static.T,
                reference.outlet_static.T);

            const double ev = relative_error(
                analytical.outlet_tangential_velocity,
                reference.outlet_tangential_velocity);

            maximum_pressure_error =
                std::max(maximum_pressure_error, ep);

            maximum_temperature_error =
                std::max(maximum_temperature_error, eT);

            maximum_velocity_error =
                std::max(maximum_velocity_error, ev);

            require(ep < 1.0e-7,
                    "Analytical pressure mismatch");

            require(eT < 1.0e-7,
                    "Analytical temperature mismatch");

            require(ev < 1.0e-7,
                    "Analytical tangential velocity mismatch");

            const double h0 =
                state.h + 0.5 * (u*u + v*v);

            require(
                std::abs(
                    analytical.stagnation_enthalpy_residual) /
                    std::max(1.0, std::abs(h0)) < 1.0e-9,
                "Analytical stagnation enthalpy defect");

            require(
                std::abs(
                    analytical.entropy_balance_residual) < 1.0e-7,
                "Analytical entropy balance defect");

            require(
                std::abs(
                    analytical.angular_momentum_balance_residual) <
                    1.0e-8,
                "Analytical angular momentum defect");
        }
    }

    std::cout
        << name
        << " max_pressure_error=" << maximum_pressure_error
        << " max_temperature_error=" << maximum_temperature_error
        << " max_velocity_error=" << maximum_velocity_error
        << '\n';

    // Comparative timing: same geometry, operating point and
    // integration resolution. This is an exploratory benchmark,
    // not a statistically controlled performance measurement.
    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = 0.02;
    options.radial_steps = 512;

    constexpr int repetitions = 20;

    volatile double sink = 0.0;

    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < repetitions; ++i) {
        const auto result =
            evaluate_frictional_vaneless_diffuser(
                model, geometry, inlet, options);

        sink = sink + result.outlet_static.p;
    }

    auto middle = std::chrono::steady_clock::now();

    for (int i = 0; i < repetitions; ++i) {
        const auto result =
            evaluate_frictional_vaneless_diffuser_analytical(
                model, model, geometry, inlet, options);

        sink = sink + result.outlet_static.p;
    }

    auto finish = std::chrono::steady_clock::now();

    const double reference_seconds =
        std::chrono::duration<double>(
            middle - start).count();

    const double analytical_seconds =
        std::chrono::duration<double>(
            finish - middle).count();

    std::cout
        << name
        << " reference_seconds=" << reference_seconds
        << " analytical_seconds=" << analytical_seconds
        << " reference_over_analytical="
        << reference_seconds / analytical_seconds
        << " sink=" << sink
        << '\n';
}

} // namespace

int main()
{
    try {
        std::cout << std::scientific
                  << std::setprecision(12);

        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        verify(ideal, "IdealGas");
        verify(pr, "PengRobinson");

        std::cout
            << "Analytical diffuser derivative comparison: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
