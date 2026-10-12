#pragma once

#include <turbocompgpu/vaneless_diffuser_branches.hpp>

#include <cmath>
#include <numbers>
#include <optional>
#include <stdexcept>

namespace turbocompgpu {

// Stage 8B.5: lossless impeller + vaneless diffuser.
//
// Station 2: impeller exit / diffuser inlet.
// Station 3: vaneless diffuser exit.
//
// A diffuser failure does not invalidate a previously
// validated impeller solution.
//
// This API does not model a collector, volute, shocks,
// diffuser losses, or compressor-system choking.

struct CompleteStageBranchResult {
    DiffuserBranch branch = DiffuserBranch::Subsonic;

    DiffuserBranchStatus status =
        DiffuserBranchStatus::BranchUnavailable;

    thermomodel::State station2_static;
    double station2_meridional_velocity = 0.0;
    double station2_tangential_velocity = 0.0;

    // Station 3 is present only when the diffuser solution
    // is converged or exactly critical.
    std::optional<DiffuserBranchPoint> station3;

    double mass_flow = 0.0;
    double impeller_specific_work = 0.0;
    double impeller_shaft_power = 0.0;

    double diffuser_critical_mass_flow = 0.0;
    double diffuser_critical_pressure = 0.0;

    // The following residuals are meaningful only when
    // station3.has_value() is true.
    double mass_residual = 0.0;
    double combined_energy_residual = 0.0;
    double diffuser_entropy_residual = 0.0;
    double diffuser_angular_momentum_residual = 0.0;

    DiffuserBranchResult diffuser;
};

inline CompleteStageBranchResult
evaluate_complete_stage_branch(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& impeller_geometry,
    const OperatingPointRequest& request,
    const OperatingPointResult& impeller,
    const VanelessDiffuserGeometry& diffuser_geometry,
    DiffuserBranch branch = DiffuserBranch::Subsonic)
{
    if (!(impeller.converged &&
          impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::invalid_argument(
            "Complete-stage impeller solution is not validated");
    }

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double r2 =
        impeller_geometry.exit.diameter / 2.0;

    if (std::abs(diffuser_geometry.inlet_radius - r2) >
        1.0e-12 * std::max(1.0, r2)) {
        throw std::invalid_argument(
            "Complete-stage impeller/diffuser radius mismatch");
    }

    const double mass = request.operating.mass_flow;

    const double area2 = effective_flow_area(
        impeller_exit_area(impeller_geometry.exit),
        impeller_geometry.exit.blockage);

    const double vm2 =
        mass / (state2.rho * area2);

    const double omega =
        2.0 * std::numbers::pi *
        request.operating.rotational_speed_rpm / 60.0;

    const double U2 =
        omega *
        impeller_geometry.impeller.exit_diameter / 2.0;

    const double beta =
        impeller_geometry.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double whirl2 =
        request.operating.slip_factor *
        (U2 - vm2 / std::tan(beta));

    const VanelessDiffuserInlet inlet{
        state2,
        vm2,
        whirl2,
        mass
    };

    CompleteStageBranchResult result;

    result.branch = branch;
    result.station2_static = state2;
    result.station2_meridional_velocity = vm2;
    result.station2_tangential_velocity = whirl2;

    result.mass_flow = mass;
    result.impeller_specific_work = impeller.specific_work;
    result.impeller_shaft_power = impeller.shaft_power;

    result.diffuser = evaluate_vaneless_diffuser_branch(
        model,
        derivatives,
        diffuser_geometry,
        inlet,
        branch);

    result.status = result.diffuser.status;

    result.diffuser_critical_mass_flow =
        result.diffuser.critical_mass_flow;

    result.diffuser_critical_pressure =
        result.diffuser.critical.pressure;

    if (result.status != DiffuserBranchStatus::Converged &&
        result.status != DiffuserBranchStatus::Critical) {
        return result;
    }

    result.station3 = result.diffuser.outlet;

    const auto state1 = evaluate(
        model,
        request.operating.inlet_temperature,
        request.operating.inlet_pressure);

    const double area1 = effective_flow_area(
        inlet_annulus_area(impeller_geometry.inlet),
        impeller_geometry.inlet.blockage);

    const double vm1 =
        mass / (state1.rho * area1);

    const double whirl1 =
        request.operating.inlet_prewhirl;

    const double h01 =
        state1.h +
        0.5 * (vm1 * vm1 + whirl1 * whirl1);

    const auto& outlet = *result.station3;

    const double h03 =
        outlet.state.h +
        0.5 * (
            outlet.meridional_velocity *
                outlet.meridional_velocity +
            outlet.tangential_velocity *
                outlet.tangential_velocity);

    result.mass_residual =
        result.diffuser.mass_residual;

    result.combined_energy_residual =
        (h03 - h01) - impeller.specific_work;

    result.diffuser_entropy_residual =
        result.diffuser.entropy_residual;

    result.diffuser_angular_momentum_residual =
        result.diffuser.angular_momentum_residual;

    return result;
}

} // namespace turbocompgpu
