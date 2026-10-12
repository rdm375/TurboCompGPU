#pragma once

#include <turbocompgpu/rotating_throat.hpp>
#include <turbocompgpu/impeller_system.hpp>
#include <turbocompgpu/thermodynamic_solver.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct InletCoupledChokingInputs {
    double inlet_temperature;
    double inlet_pressure;
    double rotational_speed_rpm;
    double inlet_prewhirl;

    double throat_radius;
    double throat_tangential_relative_velocity;
    double throat_effective_area;
};

struct InletCoupledChokingResult {
    double critical_mass_flow;
    double inlet_meridional_velocity;
    double inlet_blade_speed;
    double inlet_rothalpy;
    double effective_stagnation_temperature;
    double effective_stagnation_pressure;
    double mass_flow_residual;
    int iterations;
    RotatingThroatResult throat;
};

// Recover (T,p) from specified h and s on the vapor branch.
// Pressure is solved in logarithmic coordinates; temperature
// inversion uses TurboCompGPU's existing safeguarded solver.
inline thermomodel::State state_from_hs(
    const thermomodel::Model& model,
    double target_h,
    double target_s,
    double initial_T,
    double initial_p)
{
    const auto limits = model.metadata();

    auto evaluate_logp = [&](double logp) {
        const double p = std::exp(logp);

        const auto inversion = invert_temperature(
            model, p, target_h,
            Property::Enthalpy, initial_T);

        return inversion.state;
    };

    const double log_min = std::log(limits.pmin);
    const double log_max = std::log(limits.pmax);

    double center = std::log(initial_p);
    center = std::clamp(center, log_min, log_max);

    const auto center_state = evaluate_logp(center);
    const double center_f = center_state.s - target_s;

    if (std::abs(center_f) < 1e-11)
        return center_state;

    // At constant enthalpy, entropy decreases with pressure
    // on a stable, single-phase thermodynamic branch.
    double lo = center;
    double hi = center;
    double flo = center_f;
    double fhi = center_f;

    double step = 0.1;

    for (int i = 0; i < 100 && flo * fhi > 0.0; ++i) {
        if (center_f > 0.0) {
            const double candidate = std::min(log_max, hi + step);
            if (candidate == hi) break;
            hi = candidate;
            fhi = evaluate_logp(hi).s - target_s;
        } else {
            const double candidate = std::max(log_min, lo - step);
            if (candidate == lo) break;
            lo = candidate;
            flo = evaluate_logp(lo).s - target_s;
        }
        step *= 1.5;
    }

    if (!(flo >= 0.0 && fhi <= 0.0))
        throw std::runtime_error(
            "Effective stagnation pressure not bracketed");

    for (int i = 0; i < 120; ++i) {
        const double mid = 0.5 * (lo + hi);
        const auto state = evaluate_logp(mid);
        const double f = state.s - target_s;

        if (std::abs(f) <= 1e-10)
            return state;

        if (f > 0.0)
            lo = mid;
        else
            hi = mid;

        if (hi - lo <= 1e-13)
            break;
    }

    const auto state = evaluate_logp(0.5 * (lo + hi));

    if (std::abs(state.s - target_s) > 1e-8)
        throw std::runtime_error(
            "Effective stagnation entropy residual too large");

    return state;
}

inline InletCoupledChokingResult evaluate_inlet_coupled_choking(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& geometry,
    const InletCoupledChokingInputs& input)
{
    const auto inlet = evaluate(
        model, input.inlet_temperature, input.inlet_pressure);

    // Reuse the validated TurboCompGPU continuity geometry.
    const double inlet_area =
        effective_flow_area(
            inlet_annulus_area(geometry.inlet),
            geometry.inlet.blockage);

    if (!(std::isfinite(inlet_area) && inlet_area > 0.0))
        throw std::invalid_argument("Invalid inlet area");

    const double omega =
        2.0 * std::numbers::pi *
        input.rotational_speed_rpm / 60.0;

    const double U1 =
        omega * geometry.impeller.inlet_diameter * 0.5;

    const double Ut = omega * input.throat_radius;
    const double Wt = input.throat_tangential_relative_velocity;

    if (!(std::isfinite(U1) && std::isfinite(Ut) &&
          std::isfinite(Wt) &&
          input.throat_effective_area > 0.0))
        throw std::invalid_argument(
            "Invalid rotating throat parameters");

    auto evaluate_at_flow = [&](double mdot) {
        const double Vm1 = mdot / (inlet.rho * inlet_area);
        const double Wtheta1 = input.inlet_prewhirl - U1;

        const double I1 =
            inlet.h +
            0.5 * (Vm1 * Vm1 + Wtheta1 * Wtheta1) -
            0.5 * U1 * U1;

        const double h0eff =
            I1 + 0.5 * (Ut * Ut - Wt * Wt);

        const auto effective = state_from_hs(
            model, h0eff, inlet.s,
            inlet.T, inlet.p);

        const RotatingThroatInputs throat_input{
            effective.T,
            effective.p,
            input.throat_radius,
            input.rotational_speed_rpm,
            Wt,
            input.throat_effective_area
        };

        const auto throat = evaluate_rotating_throat(
            model, derivatives, throat_input);

        return InletCoupledChokingResult{
            mdot,
            Vm1,
            U1,
            I1,
            effective.T,
            effective.p,
            mdot - throat.critical_mass_flow,
            0,
            throat
        };
    };

    // The fixed-point residual is mdot - throat capacity.
    // Its low-flow sign must be negative for a physical
    // capacity-limited solution.
    double lo = 0.0;
    auto low = evaluate_at_flow(lo);

    if (!(low.mass_flow_residual < 0.0))
        throw std::runtime_error(
            "No positive low-flow throat capacity");

    double hi = std::max(
        1.0,
        low.throat.critical_mass_flow * 2.0);

    InletCoupledChokingResult high{};

    bool bracketed = false;

    for (int i = 0; i < 40; ++i) {
        high = evaluate_at_flow(hi);

        if (high.mass_flow_residual > 0.0) {
            bracketed = true;
            break;
        }

        hi *= 1.5;
    }

    if (!bracketed)
        throw std::runtime_error(
            "Choking mass-flow root not bracketed");

    for (int iteration = 1; iteration <= 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        auto result = evaluate_at_flow(mid);
        result.iterations = iteration;

        const double scale =
            std::max(1.0, result.critical_mass_flow);

        if (std::abs(result.mass_flow_residual) <
            1e-10 * scale)
            return result;

        if (result.mass_flow_residual < 0.0)
            lo = mid;
        else
            hi = mid;
    }

    throw std::runtime_error(
        "Inlet-coupled choking solver did not converge");
}

} // namespace turbocompgpu
