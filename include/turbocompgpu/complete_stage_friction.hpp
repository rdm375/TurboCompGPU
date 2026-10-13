#pragma once

#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/vaneless_diffuser_friction_analytical.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

// Stage 9A.1: impeller + frictional vaneless diffuser.
//
// Station 1: impeller inlet.
// Station 2: impeller exit / diffuser inlet.
// Station 3: diffuser outlet.
//
// No collector, volute, mechanical losses, or surge model.
// The diffuser has constant width and zero blockage.

struct FrictionalCompleteStageResult {
    CompleteStageStatus status = CompleteStageStatus::ImpellerFailed;

    OperatingPointResult impeller;
    std::optional<FrictionalDiffuserResult> diffuser;

    std::optional<StagnationResult> station1_stagnation;
    std::optional<StagnationResult> station3_stagnation;

    double station1_meridional_velocity = 0.0;
    double station1_tangential_velocity = 0.0;

    double station2_meridional_velocity = 0.0;
    double station2_tangential_velocity = 0.0;

    double mass_flow = 0.0;
    double fluid_specific_work = 0.0;
    double fluid_power = 0.0;

    double total_to_total_pressure_ratio = 0.0;
    double static_pressure_ratio = 0.0;

    double combined_energy_residual = 0.0;

    std::string message;

    bool complete() const noexcept {
        return status == CompleteStageStatus::Complete;
    }

    bool impeller_valid() const noexcept {
        return impeller.converged &&
               impeller.independently_validated &&
               impeller.basic_physics_valid;
    }
};

inline FrictionalCompleteStageResult
evaluate_frictional_complete_stage(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& impeller_geometry,
    const VanelessDiffuserGeometry& diffuser_geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const NonlinearSolverOptions& solver_options = {},
    const FrictionalDiffuserOptions& diffuser_options = {},
    double validation_tolerance = 1.0e-8)
{
    FrictionalCompleteStageResult result;

    result.impeller = evaluate_operating_point(
        model,
        impeller_geometry,
        request,
        scales,
        solver_options,
        validation_tolerance);

    if (!result.impeller_valid()) {
        result.status = CompleteStageStatus::ImpellerFailed;
        result.message =
            "Impeller solve or validation failed: " +
            result.impeller.message;
        return result;
    }

    try {
        const double mass = request.operating.mass_flow;

        const double r2 =
            impeller_geometry.exit.diameter / 2.0;

        if (!std::isfinite(r2) || r2 <= 0.0 ||
            !std::isfinite(diffuser_geometry.inlet_radius) ||
            std::abs(diffuser_geometry.inlet_radius - r2) >
                1.0e-12 * std::max(1.0, r2)) {
            throw std::invalid_argument(
                "Impeller/diffuser inlet radius mismatch");
        }

        if (diffuser_geometry.outlet_blockage != 0.0) {
            throw std::invalid_argument(
                "Frictional diffuser requires zero blockage");
        }

        const double area2 = effective_flow_area(
            impeller_exit_area(impeller_geometry.exit),
            impeller_geometry.exit.blockage);

        const double diffuser_inlet_area =
            2.0 * std::numbers::pi *
            diffuser_geometry.inlet_radius *
            diffuser_geometry.outlet_width;

        if (!std::isfinite(diffuser_inlet_area) ||
            !std::isfinite(area2) ||
            area2 <= 0.0 ||
            std::abs(diffuser_inlet_area - area2) >
                1.0e-10 * std::max(area2, diffuser_inlet_area)) {
            throw std::invalid_argument(
                "Impeller effective exit area does not match "
                "constant-width diffuser inlet area");
        }

        const auto state1 = evaluate(
            model,
            request.operating.inlet_temperature,
            request.operating.inlet_pressure);

        const auto state2 = evaluate(
            model,
            result.impeller.outlet_temperature,
            result.impeller.outlet_pressure);

        const double area1 = effective_flow_area(
            inlet_annulus_area(impeller_geometry.inlet),
            impeller_geometry.inlet.blockage);

        const double vm1 = mass / (state1.rho * area1);
        const double vm2 = mass / (state2.rho * area2);

        const double whirl1 =
            request.operating.inlet_prewhirl;

        const double omega =
            2.0 * std::numbers::pi *
            request.operating.rotational_speed_rpm / 60.0;

        const double U2 =
            omega * impeller_geometry.impeller.exit_diameter / 2.0;

        const double beta =
            impeller_geometry.impeller.exit_blade_angle_deg *
            std::numbers::pi / 180.0;

        const double whirl2 =
            request.operating.slip_factor *
            (U2 - vm2 / std::tan(beta));

        result.mass_flow = mass;
        result.station1_meridional_velocity = vm1;
        result.station1_tangential_velocity = whirl1;
        result.station2_meridional_velocity = vm2;
        result.station2_tangential_velocity = whirl2;

        const FlowStation station1{
            state1, {vm1, whirl1}
        };

        result.station1_stagnation =
            evaluate_stagnation_state(model, station1);

        const VanelessDiffuserInlet diffuser_inlet{
            state2, vm2, whirl2, mass
        };

        result.diffuser =
            evaluate_frictional_vaneless_diffuser_analytical(
                model,
                derivatives,
                diffuser_geometry,
                diffuser_inlet,
                diffuser_options);

        const auto& d = *result.diffuser;

        const FlowStation station3{
            d.outlet_static,
            {
                d.outlet_radial_velocity,
                d.outlet_tangential_velocity
            }
        };

        result.station3_stagnation =
            evaluate_stagnation_state(model, station3);

        const double h01 =
            result.station1_stagnation->stagnation_enthalpy;

        const double h03 =
            result.station3_stagnation->stagnation_enthalpy;

        result.fluid_specific_work = h03 - h01;
        result.fluid_power = mass * result.fluid_specific_work;

        result.combined_energy_residual =
            result.fluid_specific_work -
            result.impeller.specific_work;

        result.total_to_total_pressure_ratio =
            result.station3_stagnation->stagnation_state.p /
            result.station1_stagnation->stagnation_state.p;

        result.static_pressure_ratio =
            d.outlet_static.p / state1.p;

        if (!(std::isfinite(result.total_to_total_pressure_ratio) &&
              result.total_to_total_pressure_ratio > 0.0 &&
              std::isfinite(result.fluid_power) &&
              std::isfinite(result.combined_energy_residual))) {
            throw std::runtime_error(
                "Invalid complete-stage engineering outputs");
        }

        const double energy_scale =
            std::max(1.0, std::abs(result.impeller.specific_work));

        if (std::abs(result.combined_energy_residual) >
            1.0e-7 * energy_scale) {
            throw std::runtime_error(
                "Complete-stage energy balance failed");
        }

        result.status = CompleteStageStatus::Complete;
        result.message =
            "Frictional complete-stage calculation succeeded";
    }
    catch (const std::exception& error) {
        result.status = CompleteStageStatus::DiffuserFailed;
        result.message =
            std::string("Frictional stage failed: ") + error.what();
    }

    return result;
}

} // namespace turbocompgpu
