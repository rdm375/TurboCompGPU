#pragma once

#include <turbocompgpu/incidence_coupled_choking.hpp>
#include <turbocompgpu/throat_incidence.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Stage 6E: mass-flow-dependent prescribed transport closure.
//
// The tangential throat velocity is inferred by projecting
// the inlet relative velocity onto the throat basis.
//
// This is a frozen-relative-components transport model,
// not a predictive blade-turning or incidence-loss model.

struct PredictiveIncidenceChokingInputs {
    double inlet_temperature;
    double inlet_pressure;
    double rotational_speed_rpm;
    double inlet_prewhirl;
    ThroatGeometry throat_geometry;
};

struct PredictiveIncidenceChokingResult {
    IncidenceCoupledChokingResult choking;
    ThroatIncidenceResult incidence;

    double mass_flow_residual;
    double incidence_residual;

    int iterations;
};

inline PredictiveIncidenceChokingResult
evaluate_predictive_incidence_choking(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& geometry,
    const PredictiveIncidenceChokingInputs& input)
{
    if (!(std::isfinite(input.inlet_temperature) &&
          input.inlet_temperature > 0.0 &&
          std::isfinite(input.inlet_pressure) &&
          input.inlet_pressure > 0.0 &&
          std::isfinite(input.rotational_speed_rpm) &&
          input.rotational_speed_rpm >= 0.0 &&
          std::isfinite(input.inlet_prewhirl))) {
        throw std::invalid_argument(
            "Invalid predictive-incidence inputs");
    }

    // Validate throat geometry before root iterations.
    const auto throat_geometry =
        evaluate_throat_geometry(input.throat_geometry);

    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    if (!inlet) {
        throw std::runtime_error(
            "Predictive-incidence inlet EOS evaluation failed");
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
            incidence.transported_local_velocity.tangential;

        const IncidenceCoupledChokingInputs coupled_input{
            input.inlet_temperature,
            input.inlet_pressure,
            input.rotational_speed_rpm,
            input.inlet_prewhirl,
            input.throat_geometry,
            Wtau
        };

        const auto choking =
            evaluate_incidence_coupled_choking(
                model,
                derivatives,
                geometry,
                coupled_input);

        const double mass_residual =
            mdot - choking.coupled.critical_mass_flow;

        const double incidence_residual =
            choking.throat_local_velocity.tangential -
            incidence.transported_local_velocity.tangential;

        if (!(std::isfinite(mass_residual) &&
              std::isfinite(incidence_residual))) {
            throw std::runtime_error(
                "Nonfinite predictive-incidence residual");
        }

        return PredictiveIncidenceChokingResult{
            choking,
            incidence,
            mass_residual,
            incidence_residual,
            iteration
        };
    };

    double lo = 0.0;
    const auto low = evaluate_at_flow(lo, 0);

    if (!(low.mass_flow_residual < 0.0)) {
        throw std::runtime_error(
            "No positive low-flow choking capacity");
    }

    // Initial bracket scale from the low-flow capacity.
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
                "Nonfinite predictive-incidence bracket");
        }
    }

    if (!bracketed) {
        throw std::runtime_error(
            "Predictive-incidence mass-flow root not bracketed");
    }

    for (int iteration = 1; iteration <= 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);

        auto result =
            evaluate_at_flow(mid, iteration);

        const double scale =
            std::max(
                1.0,
                result.choking.coupled.critical_mass_flow);

        if (std::abs(result.mass_flow_residual) <
            1e-10 * scale) {

            // Reject reverse meridional flow at the sonic
            // throat. This is an admissibility check,
            // not an alternative root condition.
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
        "Predictive-incidence choking did not converge");
}

} // namespace turbocompgpu
