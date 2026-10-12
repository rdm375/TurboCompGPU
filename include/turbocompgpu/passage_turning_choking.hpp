#pragma once

#include <turbocompgpu/predictive_incidence_choking.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Stage 6F: prescribed fractional alignment with throat-normal flow.
//
// lambda = 0: blade-aligned throat (Stage 6C)
// lambda = 1: frozen-component incidence (Stage 6E)
//
// lambda is a prescribed closure parameter, not a
// calibrated aerodynamic turning or loss coefficient.

struct PassageTurningChokingInputs {
    double inlet_temperature;
    double inlet_pressure;
    double rotational_speed_rpm;
    double inlet_prewhirl;

    ThroatGeometry throat_geometry;
    double turning_fraction;
};

struct PassageTurningChokingResult {
    IncidenceCoupledChokingResult choking;
    ThroatIncidenceResult incidence;

    double turning_fraction;
    double predicted_tangential_velocity;
    double turning_residual;
    double mass_flow_residual;

    int iterations;
};

inline PassageTurningChokingResult
evaluate_passage_turning_choking(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& geometry,
    const PassageTurningChokingInputs& input)
{
    const double lambda = input.turning_fraction;

    if (!(std::isfinite(lambda) &&
          lambda >= 0.0 && lambda <= 1.0 &&
          std::isfinite(input.inlet_temperature) &&
          input.inlet_temperature > 0.0 &&
          std::isfinite(input.inlet_pressure) &&
          input.inlet_pressure > 0.0 &&
          std::isfinite(input.rotational_speed_rpm) &&
          input.rotational_speed_rpm >= 0.0 &&
          std::isfinite(input.inlet_prewhirl))) {
        throw std::invalid_argument(
            "Invalid passage-turning choking inputs");
    }

    // Validate geometry before entering the root solve.
    const auto throat_geometry =
        evaluate_throat_geometry(input.throat_geometry);

    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    if (!inlet) {
        throw std::runtime_error(
            "Passage-turning inlet EOS evaluation failed");
    }

    auto evaluate_at_flow = [&](double mdot, int iteration) {
        const auto incidence =
            evaluate_throat_incidence(
                geometry,
                {
                    inlet.value.rho,
                    mdot,
                    input.rotational_speed_rpm,
                    input.inlet_prewhirl,
                    input.throat_geometry.blade_angle_deg
                });

        const double Wtau =
            lambda *
            incidence.transported_local_velocity.tangential;

        const auto choking =
            evaluate_incidence_coupled_choking(
                model,
                derivatives,
                geometry,
                {
                    input.inlet_temperature,
                    input.inlet_pressure,
                    input.rotational_speed_rpm,
                    input.inlet_prewhirl,
                    input.throat_geometry,
                    Wtau
                });

        const double mass_residual =
            mdot - choking.coupled.critical_mass_flow;

        const double turning_residual =
            choking.throat_local_velocity.tangential - Wtau;

        if (!(std::isfinite(mass_residual) &&
              std::isfinite(turning_residual))) {
            throw std::runtime_error(
                "Nonfinite passage-turning residual");
        }

        return PassageTurningChokingResult{
            choking,
            incidence,
            lambda,
            Wtau,
            turning_residual,
            mass_residual,
            iteration
        };
    };

    double lo = 0.0;
    const auto low = evaluate_at_flow(lo, 0);

    if (!(low.mass_flow_residual < 0.0)) {
        throw std::runtime_error(
            "No positive low-flow choking capacity");
    }

    double hi = std::max(
        1.0,
        2.0 * low.choking.coupled.critical_mass_flow);

    bool bracketed = false;

    for (int i = 0; i < 40; ++i) {
        const auto high = evaluate_at_flow(hi, 0);

        if (high.mass_flow_residual > 0.0) {
            bracketed = true;
            break;
        }

        hi *= 1.5;

        if (!std::isfinite(hi)) {
            throw std::runtime_error(
                "Nonfinite passage-turning bracket");
        }
    }

    if (!bracketed) {
        throw std::runtime_error(
            "Passage-turning mass-flow root not bracketed");
    }

    for (int iteration = 1; iteration <= 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        auto result = evaluate_at_flow(mid, iteration);

        const double scale =
            std::max(
                1.0,
                result.choking.coupled.critical_mass_flow);

        if (std::abs(result.mass_flow_residual) <
            1e-10 * scale) {

            if (!(result.choking
                      .throat_relative_velocity.meridional > 0.0)) {
                throw std::runtime_error(
                    "Nonpositive sonic throat meridional flow");
            }

            if (!(throat_geometry.effective_area > 0.0)) {
                throw std::runtime_error(
                    "Invalid effective throat area");
            }

            return result;
        }

        if (result.mass_flow_residual < 0.0)
            lo = mid;
        else
            hi = mid;
    }

    throw std::runtime_error(
        "Passage-turning choking did not converge");
}

} // namespace turbocompgpu
