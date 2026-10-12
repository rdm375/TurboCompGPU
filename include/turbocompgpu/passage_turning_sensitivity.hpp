#pragma once

#include <turbocompgpu/passage_turning_choking.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Stage 6G: local sensitivity of critical mass flow to
// the prescribed passage-turning fraction.
//
// Implicit equation:
//     F(m,lambda) = m - C(lambda*g(m)) = 0
//
// C(w): Stage 6D critical mass flow at prescribed Wtau=w
// g(m): Stage 6E transported throat tangential velocity
//
// Implicit derivative:
//     dm/dlambda = C'(w)*g(m) /
//                  (1 - lambda*C'(w)*g'(m))

struct PassageTurningSensitivityResult {
    PassageTurningChokingResult operating_point;

    double capacity_derivative_w;
    double transport_derivative_m;

    double implicit_derivative_lambda;
    double finite_difference_derivative_lambda;

    double implicit_denominator;
    double finite_difference_step;
};

inline PassageTurningSensitivityResult
evaluate_passage_turning_sensitivity(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& geometry,
    const PassageTurningChokingInputs& input,
    double lambda_step = 1e-3,
    double velocity_step = 0.1)
{
    if (!(std::isfinite(lambda_step) &&
          lambda_step > 0.0 &&
          lambda_step <= 0.25 &&
          std::isfinite(velocity_step) &&
          velocity_step > 0.0)) {
        throw std::invalid_argument(
            "Invalid sensitivity difference step");
    }

    const auto operating =
        evaluate_passage_turning_choking(
            model, derivatives, geometry, input);

    const double lambda = input.turning_fraction;
    const double mdot =
        operating.choking.coupled.critical_mass_flow;

    const double w =
        operating.predicted_tangential_velocity;

    const double g =
        operating.incidence
            .transported_local_velocity.tangential;

    // Independently evaluate the linear transport slope.
    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    if (!inlet)
        throw std::runtime_error("Inlet EOS evaluation failed");

    const double area =
        effective_flow_area(
            inlet_annulus_area(geometry.inlet),
            geometry.inlet.blockage);

    const double beta =
        input.throat_geometry.blade_angle_deg *
        std::numbers::pi / 180.0;

    const double gprime =
        std::cos(beta) / (inlet.value.rho * area);

    // Stage 6D capacity derivative with respect to Wtau.
    auto capacity = [&](double wt) {
        const auto result =
            evaluate_incidence_coupled_choking(
                model, derivatives, geometry,
                {
                    input.inlet_temperature,
                    input.inlet_pressure,
                    input.rotational_speed_rpm,
                    input.inlet_prewhirl,
                    input.throat_geometry,
                    wt
                });

        return result.coupled.critical_mass_flow;
    };

    const double cprime =
        (capacity(w + velocity_step) -
         capacity(w - velocity_step)) /
        (2.0 * velocity_step);

    const double denominator =
        1.0 - lambda * cprime * gprime;

    if (!(std::isfinite(denominator) &&
          std::abs(denominator) > 1e-8)) {
        throw std::runtime_error(
            "Near-singular implicit sensitivity denominator");
    }

    const double implicit =
        cprime * g / denominator;

    auto solved_mass_flow = [&](double fraction) {
        auto changed = input;
        changed.turning_fraction = fraction;

        return evaluate_passage_turning_choking(
            model, derivatives, geometry, changed)
            .choking.coupled.critical_mass_flow;
    };

    double direct = 0.0;

    if (lambda >= lambda_step &&
        lambda <= 1.0 - lambda_step) {

        direct =
            (solved_mass_flow(lambda + lambda_step) -
             solved_mass_flow(lambda - lambda_step)) /
            (2.0 * lambda_step);

    } else if (lambda < lambda_step) {

        // Second-order forward difference, including
        // points close to the lower boundary.
        const double h =
            std::min(lambda_step, (1.0 - lambda) / 2.0);

        const double f0 = mdot;
        const double f1 = solved_mass_flow(lambda + h);
        const double f2 = solved_mass_flow(lambda + 2.0*h);

        direct = (-3.0*f0 + 4.0*f1 - f2) / (2.0*h);

    } else {

        // Second-order backward difference.
        const double h =
            std::min(lambda_step, lambda / 2.0);

        const double f0 = mdot;
        const double f1 = solved_mass_flow(lambda - h);
        const double f2 = solved_mass_flow(lambda - 2.0*h);

        direct = (3.0*f0 - 4.0*f1 + f2) / (2.0*h);
    }

    if (!(std::isfinite(implicit) &&
          std::isfinite(direct) &&
          std::isfinite(gprime))) {
        throw std::runtime_error(
            "Nonfinite passage-turning sensitivity");
    }

    return {
        operating,
        cprime,
        gprime,
        implicit,
        direct,
        denominator,
        lambda_step
    };
}

// Independently estimate the curvature at lambda=0
// using the even symmetry of the Stage 6D capacity C(w).
//
// m''(0) = C''(0) * g(m0)^2

struct PassageTurningCurvatureResult {
    double aligned_mass_flow;
    double transported_tangential_velocity;

    double capacity_second_derivative_w;
    double implicit_second_derivative_lambda;

    double direct_second_derivative_lambda;
};

inline PassageTurningCurvatureResult
evaluate_passage_turning_curvature(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& geometry,
    const PassageTurningChokingInputs& input,
    double lambda_step = 0.02,
    double velocity_step = 1.0)
{
    if (!(std::isfinite(lambda_step) &&
          lambda_step > 0.0 &&
          lambda_step <= 0.25 &&
          std::isfinite(velocity_step) &&
          velocity_step > 0.0)) {
        throw std::invalid_argument(
            "Invalid curvature difference step");
    }

    auto aligned_input = input;
    aligned_input.turning_fraction = 0.0;

    const auto aligned =
        evaluate_passage_turning_choking(
            model, derivatives, geometry, aligned_input);

    const double m0 =
        aligned.choking.coupled.critical_mass_flow;

    const double g0 =
        aligned.incidence
            .transported_local_velocity.tangential;

    auto capacity = [&](double wt) {
        return evaluate_incidence_coupled_choking(
            model, derivatives, geometry,
            {
                input.inlet_temperature,
                input.inlet_pressure,
                input.rotational_speed_rpm,
                input.inlet_prewhirl,
                input.throat_geometry,
                wt
            }).coupled.critical_mass_flow;
    };

    const double c0 = capacity(0.0);
    const double cp = capacity(velocity_step);
    const double cm = capacity(-velocity_step);

    const double csecond =
        (cp - 2.0*c0 + cm) /
        (velocity_step*velocity_step);

    const double implicit_second =
        csecond*g0*g0;

    // At the symmetric endpoint m'(0)=0, so the
    // leading direct curvature estimator is:
    //
    // 2*(m(h)-m(0))/h^2
    //
    // Richardson extrapolation removes the O(h^2)
    // truncation term for a smooth even branch.

    auto solve = [&](double lambda) {
        auto changed = input;
        changed.turning_fraction = lambda;

        return evaluate_passage_turning_choking(
            model, derivatives, geometry, changed)
            .choking.coupled.critical_mass_flow;
    };

    const double h = lambda_step;
    const double mh = solve(h);
    const double mh2 = solve(h/2.0);

    const double d2h =
        2.0*(mh - m0)/(h*h);

    const double d2h2 =
        2.0*(mh2 - m0)/((h/2.0)*(h/2.0));

    const double direct_second =
        (4.0*d2h2 - d2h)/3.0;

    if (!(std::isfinite(implicit_second) &&
          std::isfinite(direct_second))) {
        throw std::runtime_error(
            "Nonfinite passage-turning curvature");
    }

    return {
        m0,
        g0,
        csecond,
        implicit_second,
        direct_second
    };
}

} // namespace turbocompgpu
