#include <turbocompgpu/vaneless_diffuser_variable_width.hpp>
#include <turbocompgpu/vaneless_diffuser.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
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
void verify(const Model& model, const char* name)
{
    constexpr double r2 = 0.20;
    constexpr double r3 = 0.30;
    constexpr double b2 = 0.025;
    constexpr double b3 = 0.020;

    constexpr double u2 = 60.0;
    constexpr double v2 = 90.0;

    const auto state2 = evaluate(model, 340.0, 700000.0);

    const double mass =
        state2.rho *
        2.0 * std::numbers::pi * r2 * b2 * u2;

    const VanelessDiffuserInlet inlet{
        state2, u2, v2, mass
    };

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

    // Independent frictionless algebraic reference.
    // Its outlet area matches the tapered passage at r3.
    //
    // This reference uses the existing pressure-root solver,
    // not the frictional diffuser's radial ODE.

    const VanelessDiffuserGeometry reference_geometry{
        r2, r3, b3, 0.0
    };

    const auto reference =
        evaluate_vaneless_diffuser(
            model, reference_geometry, inlet);

    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = 0.0;
    options.radial_steps = 512;

    const auto frictionless =
        evaluate_variable_width_vaneless_diffuser(
            model, model, tapered, inlet, options);

    const double ep = relative_error(
        frictionless.outlet_static.p,
        reference.outlet_static.p);

    const double eT = relative_error(
        frictionless.outlet_static.T,
        reference.outlet_static.T);

    const double eu = relative_error(
        frictionless.outlet_radial_velocity,
        reference.outlet_meridional_velocity);

    const double ev = relative_error(
        frictionless.outlet_tangential_velocity,
        reference.outlet_tangential_velocity);

    std::cout
        << name
        << " frictionless_reference"
        << " pressure_error=" << ep
        << " temperature_error=" << eT
        << " radial_velocity_error=" << eu
        << " tangential_velocity_error=" << ev
        << '\n';

    require(ep < 1.0e-8,
            "Independent frictionless pressure mismatch");

    require(eT < 1.0e-8,
            "Independent frictionless temperature mismatch");

    require(eu < 1.0e-8,
            "Independent frictionless radial velocity mismatch");

    require(ev < 1.0e-8,
            "Independent frictionless tangential velocity mismatch");

    // Frictional tapered-passage convergence.
    // Compare successively refined solutions.
    //
    // For classical RK4, the asymptotic ratio
    // |y_N-y_2N| / |y_2N-y_4N| approaches 16.

    options.darcy_friction_factor = 0.02;

    constexpr std::array<int, 5> meshes{
        8, 16, 32, 64, 128
    };

    std::array<double, meshes.size()> pressures{};
    std::array<double, meshes.size()> temperatures{};

    for (std::size_t i = 0; i < meshes.size(); ++i) {
        options.radial_steps = meshes[i];

        const auto result =
            evaluate_variable_width_vaneless_diffuser(
                model, model, tapered, inlet, options);

        pressures[i] = result.outlet_static.p;
        temperatures[i] = result.outlet_static.T;

        std::cout
            << name
            << " mesh=" << meshes[i]
            << " p3=" << result.outlet_static.p
            << " T3=" << result.outlet_static.T
            << " energy_residual="
            << result.stagnation_enthalpy_residual
            << " entropy_residual="
            << result.entropy_balance_residual
            << " angular_residual="
            << result.angular_momentum_balance_residual
            << '\n';
    }

    auto observed_order = [](double a, double b, double c) {
        const double numerator = std::abs(a - b);
        const double denominator = std::abs(b - c);

        if (!(numerator > 0.0 && denominator > 0.0))
            throw std::runtime_error(
                "Convergence differences are zero");

        return std::log2(numerator / denominator);
    };

    // Use the coarser triplets to reduce sensitivity
    // to floating-point roundoff at high resolution.

    const double qp = observed_order(
        pressures[0], pressures[1], pressures[2]);

    const double qT = observed_order(
        temperatures[0], temperatures[1], temperatures[2]);

    std::cout
        << name
        << " tapered_RK4"
        << " pressure_order=" << qp
        << " temperature_order=" << qT
        << '\n';

    require(qp > 3.5 && qp < 4.5,
            "Pressure convergence is not fourth order");

    require(qT > 3.5 && qT < 4.5,
            "Temperature convergence is not fourth order");
}

} // namespace

int main()
{
    try {
        std::cout
            << std::scientific
            << std::setprecision(12);

        const thermomodel::IdealGas ideal;

        verify(ideal, "IdealGas");

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        verify(pr, "PengRobinson");

        std::cout
            << "Variable-width independent verification: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
