#include <turbocompgpu/vaneless_diffuser_friction.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void require(bool condition, const char* message)
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
void test_model(const Model& model, const char* name)
{
    const VanelessDiffuserGeometry geometry{
        0.20, 0.24, 0.025, 0.0
    };

    const auto state2 = evaluate(model, 340.0, 700000.0);

    const double u2 = 60.0;
    const double v2 = 90.0;

    const double area2 =
        2.0 * std::numbers::pi *
        geometry.inlet_radius * geometry.outlet_width;

    const VanelessDiffuserInlet inlet{
        state2,
        u2,
        v2,
        state2.rho * area2 * u2
    };

    // Zero-friction comparison with the existing lossless solver.
    const auto reference =
        evaluate_vaneless_diffuser(model, geometry, inlet);

    FrictionalDiffuserOptions options;
    options.darcy_friction_factor = 0.0;
    options.radial_steps = 512;

    const auto zero =
        evaluate_frictional_vaneless_diffuser(
            model, geometry, inlet, options);

    const double pressure_error =
        relative_error(
            zero.outlet_static.p,
            reference.outlet_static.p);

    const double velocity_error =
        relative_error(
            zero.outlet_radial_velocity,
            reference.outlet_meridional_velocity);

    require(pressure_error < 1.0e-5,
            "Zero-friction pressure mismatch");

    require(velocity_error < 1.0e-5,
            "Zero-friction velocity mismatch");

    // Frictional calculation.
    options.darcy_friction_factor = 0.02;

    const auto frictional =
        evaluate_frictional_vaneless_diffuser(
            model, geometry, inlet, options);

    require(frictional.entropy_change > 0.0,
            "Friction did not increase entropy");

    require(frictional.integrated_entropy_production > 0.0,
            "Integrated entropy production is not positive");

    require(frictional.angular_momentum_change < 0.0,
            "Friction did not reduce angular momentum");

    require(
        relative_error(
            frictional.stagnation_enthalpy_residual,
            0.0) < 1.0e-5,
        "Stagnation enthalpy conservation failed");

    require(
        std::abs(frictional.entropy_balance_residual) < 1.0e-3,
        "Entropy balance failed");

    require(
        std::abs(
            frictional.angular_momentum_balance_residual) < 1.0e-3,
        "Angular momentum balance failed");

    // Radial refinement.
    options.radial_steps = 1024;

    const auto refined =
        evaluate_frictional_vaneless_diffuser(
            model, geometry, inlet, options);

    require(
        relative_error(
            frictional.outlet_static.p,
            refined.outlet_static.p) < 1.0e-6,
        "Radial refinement failed");

    std::cout
        << name
        << " zero_friction_pressure_error=" << pressure_error
        << " zero_friction_velocity_error=" << velocity_error
        << " entropy_change=" << frictional.entropy_change
        << " entropy_balance=" <<
            frictional.entropy_balance_residual
        << " angular_momentum_balance=" <<
            frictional.angular_momentum_balance_residual
        << " energy_residual=" <<
            frictional.stagnation_enthalpy_residual
        << '\n';
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        test_model(ideal, "IdealGas");
        test_model(pr, "PengRobinson");

        std::cout
            << "Frictional diffuser initial verification: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
