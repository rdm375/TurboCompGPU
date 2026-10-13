#include <turbocompgpu/complete_stage_friction.hpp>
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
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

double relative_error(double actual, double expected)
{
    return std::abs(actual - expected) /
        std::max({1.0, std::abs(actual), std::abs(expected)});
}

void check_close(
    const std::string& name,
    double actual,
    double expected,
    double tolerance)
{
    const double error = relative_error(actual, expected);

    if (!std::isfinite(error) || error > tolerance) {
        std::cerr << name
                  << " actual=" << actual
                  << " expected=" << expected
                  << " relative_error=" << error
                  << " tolerance=" << tolerance << '\n';
        throw std::runtime_error(name);
    }
}

struct Fixture {
    ImpellerSystemGeometry impeller_geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    NonlinearResidualScales scales{
        50000.0, 100.0, 500000.0
    };

    OperatingPointRequest request{
        {
            300.0,
            500000.0,
            5.0,
            12000.0,
            0.90,
            20.0,
            20.0
        },
        {
            340.0,
            std::log(700000.0 / scales.pressure)
        }
    };

    NonlinearSolverOptions solver_options;
    VanelessDiffuserGeometry diffuser_geometry;

    Fixture()
    {
        solver_options.residual_tolerance = 1.0e-13;

        const double r2 =
            impeller_geometry.exit.diameter / 2.0;

        const double area2 = effective_flow_area(
            impeller_exit_area(impeller_geometry.exit),
            impeller_geometry.exit.blockage);

        const double width =
            area2 / (2.0 * std::numbers::pi * r2);

        diffuser_geometry = {
            r2,
            1.5 * r2,
            width,
            0.0
        };
    }
};

template<class Model>
FrictionalCompleteStageResult solve(
    const Model& model,
    const Fixture& fixture,
    double friction,
    int radial_steps)
{
    FrictionalDiffuserOptions diffuser_options;
    diffuser_options.darcy_friction_factor = friction;
    diffuser_options.radial_steps = radial_steps;

    return evaluate_frictional_complete_stage(
        model,
        model,
        fixture.impeller_geometry,
        fixture.diffuser_geometry,
        fixture.request,
        fixture.scales,
        fixture.solver_options,
        diffuser_options);
}

template<class Model>
void verify_frictionless_limit(
    const Model& model,
    const char* name,
    const Fixture& fixture)
{
    const auto stage = solve(model, fixture, 0.0, 1024);

    require(stage.complete(),
            std::string(name) + ": frictionless stage failed: " +
            stage.message);

    const VanelessDiffuserInlet inlet{
        stage.diffuser->inlet_static,
        stage.station2_meridional_velocity,
        stage.station2_tangential_velocity,
        stage.mass_flow
    };

    // Independent algebraic lossless reference.
    const auto reference = evaluate_vaneless_diffuser(
        model, fixture.diffuser_geometry, inlet);

    const auto& d = *stage.diffuser;

    check_close(
        "frictionless pressure",
        d.outlet_static.p,
        reference.outlet_static.p,
        1.0e-8);

    check_close(
        "frictionless temperature",
        d.outlet_static.T,
        reference.outlet_static.T,
        1.0e-8);

    check_close(
        "frictionless radial velocity",
        d.outlet_radial_velocity,
        reference.outlet_meridional_velocity,
        1.0e-8);

    check_close(
        "frictionless tangential velocity",
        d.outlet_tangential_velocity,
        reference.outlet_tangential_velocity,
        1.0e-10);

    check_close(
        "frictionless entropy",
        d.entropy_change,
        0.0,
        1.0e-7);

    check_close(
        "frictionless angular momentum",
        d.angular_momentum_change,
        0.0,
        1.0e-8);

    check_close(
        "frictionless stagnation pressure",
        stage.station3_stagnation->stagnation_state.p,
        reference.outlet_stagnation_pressure,
        1.0e-8);

    std::cout
        << name
        << " frictionless_limit=PASS"
        << " pressure_error="
        << relative_error(
            d.outlet_static.p,
            reference.outlet_static.p)
        << '\n';
}

template<class Model>
void verify_grid_convergence(
    const Model& model,
    const char* name,
    const Fixture& fixture)
{
    constexpr std::array<int, 6> meshes{
        16, 32, 64, 128, 256, 512
    };

    std::array<double, meshes.size()> pressures{};

    std::cout << "\n" << name << " GRID CONVERGENCE\n";

    for (std::size_t i = 0; i < meshes.size(); ++i) {
        const auto stage = solve(
            model, fixture, 0.01, meshes[i]);

        require(stage.complete(),
                std::string(name) +
                ": grid solve failed at N=" +
                std::to_string(meshes[i]) +
                ": " + stage.message);

        pressures[i] = stage.diffuser->outlet_static.p;

        std::cout
            << "N=" << meshes[i]
            << " p3=" << pressures[i]
            << " energy_residual="
            << stage.combined_energy_residual
            << '\n';
    }

    // Consecutive differences provide an observed-order
    // estimate without requiring an exact solution.
    //
    // At very fine meshes, floating-point noise can dominate.
    // We therefore test fourth-order behavior on the coarser
    // resolved meshes, and separately test final stabilization.

    const double d16_32 =
        std::abs(pressures[0] - pressures[1]);

    const double d32_64 =
        std::abs(pressures[1] - pressures[2]);

    const double d64_128 =
        std::abs(pressures[2] - pressures[3]);

    require(d16_32 > 0.0 && d32_64 > 0.0 &&
            d64_128 > 0.0,
            "Grid differences are unresolved");

    const double order1 = std::log2(d16_32 / d32_64);
    const double order2 = std::log2(d32_64 / d64_128);

    std::cout
        << "observed_order_1=" << order1
        << " observed_order_2=" << order2
        << '\n';

    require(
        order1 > 3.5 && order1 < 4.5 &&
        order2 > 3.5 && order2 < 4.5,
        std::string(name) +
        ": expected fourth-order RK4 convergence");

    require(
        relative_error(pressures[4], pressures[5]) <
            1.0e-8,
        std::string(name) +
        ": pressure failed to stabilize");
}

template<class Model>
void verify_friction_sensitivity(
    const Model& model,
    const char* name,
    const Fixture& fixture)
{
    constexpr std::array<double, 5> factors{
        0.0, 0.0025, 0.005, 0.01, 0.02
    };

    double previous_entropy = -1.0;
    double previous_p03 = 0.0;
    double previous_L3 = 0.0;
    bool first = true;

    std::cout << "\n" << name << " FRICTION SWEEP\n";

    for (double friction : factors) {
        const auto stage = solve(
            model, fixture, friction, 512);

        require(stage.complete(),
                std::string(name) +
                ": friction sweep failed at f_D=" +
                std::to_string(friction) +
                ": " + stage.message);

        const auto& d = *stage.diffuser;

        const double p03 =
            stage.station3_stagnation->stagnation_state.p;

        const double L3 =
            fixture.diffuser_geometry.outlet_radius *
            d.outlet_tangential_velocity;

        require(d.entropy_change >= -1.0e-7,
                "Negative entropy generation");

        require(d.integrated_entropy_production >= -1.0e-12,
                "Negative integrated entropy production");

        check_close(
            "entropy integral",
            d.entropy_change,
            d.integrated_entropy_production,
            1.0e-7);

        check_close(
            "angular momentum integral",
            d.angular_momentum_change,
            d.integrated_wall_torque_per_mass,
            1.0e-7);

        check_close(
            "stage energy",
            stage.combined_energy_residual,
            0.0,
            1.0e-7);

        if (!first) {
            require(
                d.entropy_change >= previous_entropy - 1.0e-7,
                "Entropy not monotonic over friction sweep");

            require(
                p03 <= previous_p03 +
                    1.0e-8 * std::max(1.0, previous_p03),
                "Stagnation pressure increased with friction");

            require(
                L3 <= previous_L3 + 1.0e-8,
                "Angular momentum increased with friction");
        }

        std::cout
            << "f_D=" << friction
            << " p3=" << d.outlet_static.p
            << " p03=" << p03
            << " ds=" << d.entropy_change
            << " L3=" << L3
            << '\n';

        previous_entropy = d.entropy_change;
        previous_p03 = p03;
        previous_L3 = L3;
        first = false;
    }
}

template<class Model>
void run_model(const Model& model, const char* name)
{
    const Fixture fixture;

    verify_frictionless_limit(model, name, fixture);
    verify_grid_convergence(model, name, fixture);
    verify_friction_sensitivity(model, name, fixture);
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
            << "\nStage 9A.2 independent verification: PASS\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
