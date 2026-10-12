#pragma once

#include <turbocompgpu/passage_turning_choking.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace turbocompgpu {

// Stage 6H: diagnostic continuation of the choking boundary.
//
// The previous solution is retained as a predictor, but
// the validated Stage 6F bracketed solver remains authoritative.
//
// A small implicit denominator signals loss of local
// invertibility of F(m,lambda)=0; it does not establish
// a physical compressor-map fold.

enum class ChokingSweepParameter {
    rotational_speed_rpm,
    inlet_prewhirl,
    throat_blade_angle_deg,
    turning_fraction
};

struct ChokingContinuationPoint {
    double parameter_value;
    double predicted_mass_flow;
    double critical_mass_flow;
    double predictor_error;

    double tangential_velocity;
    double transported_tangential_velocity;

    double implicit_denominator;
    double capacity_derivative_w;
    double transport_derivative_m;

    double mass_flow_residual;
    double rothalpy_residual;
    double continuity_residual;
    double sonic_residual;

    double throat_meridional_velocity;
    int iterations;

    bool near_singular;
};

struct ChokingContinuationOptions {
    double velocity_difference_step = 0.1;
    double singularity_threshold = 0.05;
};

inline PassageTurningChokingInputs
with_choking_sweep_parameter(
    PassageTurningChokingInputs input,
    ChokingSweepParameter parameter,
    double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("Nonfinite sweep parameter");

    switch (parameter) {
    case ChokingSweepParameter::rotational_speed_rpm:
        input.rotational_speed_rpm = value;
        break;
    case ChokingSweepParameter::inlet_prewhirl:
        input.inlet_prewhirl = value;
        break;
    case ChokingSweepParameter::throat_blade_angle_deg:
        input.throat_geometry.blade_angle_deg = value;
        break;
    case ChokingSweepParameter::turning_fraction:
        input.turning_fraction = value;
        break;
    }

    return input;
}

template<class Model>
std::vector<ChokingContinuationPoint>
evaluate_choking_continuation(
    const Model& model,
    const ImpellerSystemGeometry& geometry,
    const PassageTurningChokingInputs& baseline,
    ChokingSweepParameter parameter,
    const std::vector<double>& values,
    ChokingContinuationOptions options = {})
{
    if (!(std::isfinite(options.velocity_difference_step) &&
          options.velocity_difference_step > 0.0 &&
          std::isfinite(options.singularity_threshold) &&
          options.singularity_threshold > 0.0)) {
        throw std::invalid_argument(
            "Invalid choking-continuation options");
    }

    std::vector<ChokingContinuationPoint> points;
    points.reserve(values.size());

    for (double value : values) {
        const auto input =
            with_choking_sweep_parameter(
                baseline, parameter, value);

        const auto result =
            evaluate_passage_turning_choking(
                model, model, geometry, input);

        const double mdot =
            result.choking.coupled.critical_mass_flow;

        const double w =
            result.predicted_tangential_velocity;

        const double lambda = input.turning_fraction;

        const double h = options.velocity_difference_step;

        auto capacity = [&](double wt) {
            return evaluate_incidence_coupled_choking(
                model, model, geometry,
                {
                    input.inlet_temperature,
                    input.inlet_pressure,
                    input.rotational_speed_rpm,
                    input.inlet_prewhirl,
                    input.throat_geometry,
                    wt
                }).coupled.critical_mass_flow;
        };

        const double cprime =
            (capacity(w+h)-capacity(w-h))/(2.0*h);

        const auto inlet = model.evaluate_Tp(
            input.inlet_temperature,
            input.inlet_pressure);

        if (!inlet)
            throw std::runtime_error(
                "Continuation inlet EOS evaluation failed");

        const double area =
            effective_flow_area(
                inlet_annulus_area(geometry.inlet),
                geometry.inlet.blockage);

        const double beta =
            input.throat_geometry.blade_angle_deg *
            std::numbers::pi/180.0;

        const double gprime =
            std::cos(beta)/(inlet.value.rho*area);

        const double denominator =
            1.0-lambda*cprime*gprime;

        const double predictor =
            points.empty()
                ? mdot
                : points.back().critical_mass_flow;

        const double predictor_error =
            mdot-predictor;

        if (!(std::isfinite(mdot) &&
              std::isfinite(cprime) &&
              std::isfinite(gprime) &&
              std::isfinite(denominator) &&
              std::isfinite(predictor_error))) {
            throw std::runtime_error(
                "Nonfinite continuation diagnostic");
        }

        points.push_back({
            value,
            predictor,
            mdot,
            predictor_error,
            w,
            result.incidence
                .transported_local_velocity.tangential,
            denominator,
            cprime,
            gprime,
            result.mass_flow_residual,
            result.choking.rothalpy_residual,
            result.choking.continuity_residual,
            result.choking.sonic_residual,
            result.choking.throat_relative_velocity.meridional,
            result.iterations,
            std::abs(denominator) <
                options.singularity_threshold
        });
    }

    return points;
}

} // namespace turbocompgpu
