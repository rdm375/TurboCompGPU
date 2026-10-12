#include <turbocompgpu/vaneless_diffuser.hpp>
#include <turbocompgpu/complete_stage_branches.hpp>

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

void check_close(
    const char* name,
    double actual,
    double expected,
    double relative_tolerance = 1.0e-8)
{
    const double scale = std::max({
        1.0, std::abs(actual), std::abs(expected)
    });

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) >
            relative_tolerance * scale) {
        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected << '\n';

        throw std::runtime_error(name);
    }
}

void check(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template <typename ModelType>
void run_case(
    const ModelType& model,
    const char* model_name)
{
    const ImpellerSystemGeometry geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearResidualScales scales{
        50000.0, 100.0, 500000.0
    };

    const OperatingPointRequest request{
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

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    // Step 1: solve the actual impeller operating point.
    const auto impeller = evaluate_operating_point(
        model, geometry, request, scales, options);

    check(
        impeller.converged &&
        impeller.independently_validated &&
        impeller.basic_physics_valid,
        "Impeller operating point failed validation"
    );

    const auto state1 = evaluate(
        model,
        request.operating.inlet_temperature,
        request.operating.inlet_pressure);

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double m = request.operating.mass_flow;

    const double area1 = effective_flow_area(
        inlet_annulus_area(geometry.inlet),
        geometry.inlet.blockage);

    const double area2 = effective_flow_area(
        impeller_exit_area(geometry.exit),
        geometry.exit.blockage);

    const double vm1 = m / (state1.rho * area1);
    const double vm2 = m / (state2.rho * area2);

    const double omega =
        2.0 * std::numbers::pi *
        request.operating.rotational_speed_rpm / 60.0;

    const double r1 =
        geometry.impeller.inlet_diameter / 2.0;

    const double r2 =
        geometry.impeller.exit_diameter / 2.0;

    const double beta =
        geometry.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double whirl1 =
        request.operating.inlet_prewhirl;

    const double whirl2 =
        request.operating.slip_factor *
        (omega * r2 - vm2 / std::tan(beta));

    // Independently reconstruct the impeller energy balance.
    const double h01 =
        state1.h + 0.5 * (vm1 * vm1 + whirl1 * whirl1);

    const double h02 =
        state2.h + 0.5 * (vm2 * vm2 + whirl2 * whirl2);

    const double euler_work =
        omega * (r2 * whirl2 - r1 * whirl1);

    check_close(
        "Impeller Euler work",
        euler_work,
        impeller.specific_work);

    check_close(
        "Impeller energy conservation",
        h02 - h01,
        euler_work);

    check_close(
        "Impeller power",
        m * euler_work,
        impeller.shaft_power);

    // Step 2: construct a feasible, fixed diffuser geometry.
    //
    // The geometry is manufactured from a known isentropic
    // state, then passed to the predictive diffuser solver.
    // This avoids guessing a width that may admit no solution.
    //
    // The resulting diffuser still solves for pressure:
    // its outlet pressure is NOT supplied to the solver.

    const FlowStation impeller_exit_station{
        state2, {vm2, whirl2}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(
            model, impeller_exit_station);

    const double p02 =
        stagnation2.stagnation_state.p;

    check(p02 > state2.p,
          "Invalid impeller-exit stagnation pressure");

    const double r3 = 1.5 * r2;

    const double whirl3 =
        r2 * whirl2 / r3;

    // Choose a reference pressure near the high-pressure
    // recovery branch, but below stagnation pressure.
    const double p3_reference =
        state2.p + 0.10 * (p02 - state2.p);

    const auto reference = invert_temperature(
        model,
        p3_reference,
        state2.s,
        Property::Entropy,
        state2.T);

    const double kinetic3 =
        2.0 * (h02 - reference.state.h) -
        whirl3 * whirl3;

    check(kinetic3 > 0.0,
          "Manufactured diffuser kinetic energy invalid");

    const double vm3_reference =
        std::sqrt(kinetic3);

    const double area3 =
        m / (reference.state.rho * vm3_reference);

    const double width3 =
        area3 / (2.0 * std::numbers::pi * r3);

    const VanelessDiffuserGeometry diffuser_geometry{
        r2,
        r3,
        width3,
        0.0
    };

    // Step 3: exercise the actual impeller-to-diffuser adapter.
    const auto diffuser =
        evaluate_vaneless_diffuser_from_operating_point(
            model,
            geometry,
            request,
            impeller,
            diffuser_geometry);

    // Step 4: verify the interface velocity reconstruction.
    check_close(
        "Impeller/diffuser meridional interface",
        diffuser.inlet_meridional_velocity,
        vm2);

    check_close(
        "Impeller/diffuser whirl interface",
        diffuser.inlet_tangential_velocity,
        whirl2);

    check_close(
        "Impeller/diffuser pressure interface",
        diffuser.inlet_static.p,
        state2.p);

    check_close(
        "Impeller/diffuser temperature interface",
        diffuser.inlet_static.T,
        state2.T);

    // Step 5: verify diffuser conservation independently.
    const double h03 =
        diffuser.outlet_static.h +
        0.5 * (
            diffuser.outlet_meridional_velocity *
                diffuser.outlet_meridional_velocity +
            diffuser.outlet_tangential_velocity *
                diffuser.outlet_tangential_velocity
        );

    const double reconstructed_mass =
        diffuser.outlet_static.rho *
        diffuser.outlet_effective_area *
        diffuser.outlet_meridional_velocity;

    check_close(
        "Diffuser mass conservation",
        reconstructed_mass,
        m);

    check_close(
        "Diffuser energy conservation",
        h03,
        h02);

    check_close(
        "Diffuser entropy conservation",
        diffuser.outlet_static.s,
        state2.s);

    check_close(
        "Diffuser angular momentum conservation",
        r3 * diffuser.outlet_tangential_velocity,
        r2 * whirl2);

    check_close(
        "Diffuser stagnation pressure",
        diffuser.outlet_stagnation_pressure,
        diffuser.inlet_stagnation_pressure);

    // Step 6: verify the combined impeller-diffuser balance.
    check_close(
        "Combined energy conservation",
        h03 - h01,
        euler_work);

    check_close(
        "Combined fluid power",
        m * (h03 - h01),
        impeller.shaft_power);

    check(
        diffuser.outlet_static.p > state2.p,
        "Diffuser did not recover static pressure");

    check(
        diffuser.admissible_roots >= 1,
        "Diffuser did not report an admissible root");


    // Stage 8B.5: branch-aware complete-stage integration.

    const auto stage = evaluate_complete_stage_branch(
        model, model,
        geometry, request, impeller,
        diffuser_geometry);

    check(
        stage.status == DiffuserBranchStatus::Converged,
        "Complete-stage subsonic solve failed");

    check(
        stage.station3.has_value(),
        "Complete-stage subsonic outlet missing");

    check(
        stage.branch == DiffuserBranch::Subsonic,
        "Default complete-stage branch is not subsonic");

    check_close(
        "Station 2 pressure",
        stage.station2_static.p,
        state2.p);

    check_close(
        "Station 2 meridional velocity",
        stage.station2_meridional_velocity,
        vm2);

    check_close(
        "Station 2 tangential velocity",
        stage.station2_tangential_velocity,
        whirl2);

    check_close(
        "Station 3 pressure",
        stage.station3->pressure,
        diffuser.outlet_static.p);

    check_close(
        "Station 3 temperature",
        stage.station3->state.T,
        diffuser.outlet_static.T);

    check(
        stage.station3->meridional_mach < 1.0,
        "Subsonic Mach classification failed");

    check_close(
        "Complete-stage mass conservation",
        stage.mass_residual,
        0.0,
        1.0e-8);

    check_close(
        "Complete-stage energy conservation",
        stage.combined_energy_residual,
        0.0,
        1.0e-7);

    check_close(
        "Complete-stage entropy conservation",
        stage.diffuser_entropy_residual,
        0.0,
        1.0e-7);

    check_close(
        "Complete-stage angular momentum conservation",
        stage.diffuser_angular_momentum_residual,
        0.0,
        1.0e-7);

    // Explicitly request the supersonic branch.

    const auto supersonic_stage =
        evaluate_complete_stage_branch(
            model, model,
            geometry, request, impeller,
            diffuser_geometry,
            DiffuserBranch::Supersonic);

    check(
        supersonic_stage.status ==
            DiffuserBranchStatus::Converged,
        "Complete-stage supersonic solve failed");

    check(
        supersonic_stage.station3.has_value(),
        "Complete-stage supersonic outlet missing");

    check(
        supersonic_stage.station3->meridional_mach > 1.0,
        "Supersonic Mach classification failed");

    check(
        supersonic_stage.station3->pressure <
            stage.diffuser_critical_pressure,
        "Supersonic pressure ordering failed");

    check(
        stage.diffuser_critical_pressure <
            stage.station3->pressure,
        "Subsonic pressure ordering failed");

    check_close(
        "Complete-stage critical capacity",
        stage.diffuser_critical_mass_flow,
        supersonic_stage.diffuser_critical_mass_flow,
        1.0e-10);

    check_close(
        "Complete-stage supersonic energy conservation",
        supersonic_stage.combined_energy_residual,
        0.0,
        1.0e-7);

    std::cout
        << model_name
        << " complete stage:"
        << " p2=" << stage.station2_static.p
        << " p3_sub=" << stage.station3->pressure
        << " p3_super=" << supersonic_stage.station3->pressure
        << " m_critical=" << stage.diffuser_critical_mass_flow
        << " energy_residual=" << stage.combined_energy_residual
        << '\n';

    std::cout
        << model_name
        << ": p2=" << state2.p
        << " p3=" << diffuser.outlet_static.p
        << " roots=" << diffuser.admissible_roots
        << " work=" << euler_work
        << " energy_residual=" << diffuser.energy_residual
        << '\n';

    std::cout
        << "PASS: " << model_name
        << " end-to-end impeller/diffuser integration\n";
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas ideal;
        run_case(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters params;
        params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(caloric, params);
        run_case(pr, "PengRobinson");

        std::cout
            << "PASS: Stage 8B.3 integration\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
