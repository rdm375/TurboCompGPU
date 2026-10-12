#pragma once

#include <turbocompgpu/complete_stage.hpp>
#include <turbocompgpu/continuity.hpp>
#include <turbocompgpu/operating_point.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace turbocompgpu {

// Stage 8B.2: fixed-geometry, adiabatic, lossless vaneless diffuser.
//
// Conservation equations:
//   m = rho * A * Vm
//   h + (Vm^2 + Vtheta^2)/2 = h02
//   s = s2
//   r * Vtheta = constant
//
// Pressure roots are located by scanning the admissible
// pressure interval and refining sign-changing brackets.
//
// No empirical diffuser losses or choking correlations are used.

struct VanelessDiffuserGeometry {
    double inlet_radius;
    double outlet_radius;
    double outlet_width;
    double outlet_blockage;
};

struct VanelessDiffuserInlet {
    thermomodel::State static_state;
    double meridional_velocity;
    double tangential_velocity;
    double mass_flow;
};

struct VanelessDiffuserResult {
    thermomodel::State inlet_static;
    thermomodel::State outlet_static;

    double inlet_meridional_velocity;
    double inlet_tangential_velocity;
    double outlet_meridional_velocity;
    double outlet_tangential_velocity;

    double outlet_effective_area;
    double outlet_absolute_speed;

    double inlet_stagnation_pressure;
    double outlet_stagnation_pressure;

    double mass_residual;
    double energy_residual;
    double entropy_residual;
    double angular_momentum_residual;

    int admissible_roots;
    double selected_pressure;
};

inline VanelessDiffuserResult evaluate_vaneless_diffuser(
    const thermomodel::Model& model,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet)
{
    if (!(std::isfinite(geometry.inlet_radius) &&
          std::isfinite(geometry.outlet_radius) &&
          std::isfinite(geometry.outlet_width) &&
          std::isfinite(geometry.outlet_blockage) &&
          geometry.inlet_radius > 0.0 &&
          geometry.outlet_radius > geometry.inlet_radius &&
          geometry.outlet_width > 0.0 &&
          geometry.outlet_blockage >= 0.0 &&
          geometry.outlet_blockage < 1.0 &&
          std::isfinite(inlet.meridional_velocity) &&
          inlet.meridional_velocity > 0.0 &&
          std::isfinite(inlet.tangential_velocity) &&
          std::isfinite(inlet.mass_flow) &&
          inlet.mass_flow > 0.0)) {
        throw std::invalid_argument(
            "Invalid vaneless diffuser geometry or inlet");
    }

    const auto state2 = evaluate(
        model,
        inlet.static_state.T,
        inlet.static_state.p);

    const double area3 =
        2.0 * std::numbers::pi *
        geometry.outlet_radius *
        geometry.outlet_width *
        (1.0 - geometry.outlet_blockage);

    const double whirl3 =
        geometry.inlet_radius /
        geometry.outlet_radius *
        inlet.tangential_velocity;

    const FlowStation station2{
        state2,
        {inlet.meridional_velocity,
         inlet.tangential_velocity}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, station2);

    const double h02 =
        stagnation2.stagnation_enthalpy;

    const double s2 = state2.s;
    const double p2 = state2.p;
    const double p02 = stagnation2.stagnation_state.p;

    if (!(p02 > p2)) {
        throw std::runtime_error(
            "No admissible diffuser pressure interval");
    }

    struct Trial {
        thermomodel::State state;
        double meridional_velocity;
        double residual;
    };

    auto trial_at = [&](double pressure) -> Trial {
        const auto inversion = invert_temperature(
            model,
            pressure,
            s2,
            Property::Entropy,
            state2.T);

        const auto& state = inversion.state;

        const double vm =
            inlet.mass_flow / (state.rho * area3);

        const double residual =
            state.h +
            0.5 * (vm * vm + whirl3 * whirl3) -
            h02;

        if (!std::isfinite(residual)) {
            throw std::runtime_error(
                "Non-finite diffuser energy residual");
        }

        return {state, vm, residual};
    };

    // The scan is deterministic and searches the full
    // static-to-stagnation pressure interval.
    //
    // A root is not accepted solely because an iteration
    // converged: it must satisfy the governing equations.

    constexpr int intervals = 512;
    constexpr int iterations = 100;

    std::vector<double> roots;

    double previous_pressure = p2;
    auto previous = trial_at(previous_pressure);

    const double energy_scale =
        std::max(1.0, std::abs(h02));

    for (int i = 1; i <= intervals; ++i) {
        const double fraction =
            static_cast<double>(i) / intervals;

        const double pressure =
            p2 + fraction * (p02 - p2);

        const auto current = trial_at(pressure);

        if ((previous.residual < 0.0) !=
            (current.residual < 0.0)) {

            double lo = previous_pressure;
            double hi = pressure;
            double flo = previous.residual;

            for (int j = 0; j < iterations; ++j) {
                const double mid =
                    lo + 0.5 * (hi - lo);

                if (mid == lo || mid == hi)
                    break;

                const double fm = trial_at(mid).residual;

                if ((fm < 0.0) == (flo < 0.0)) {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
            }

            const double root =
                lo + 0.5 * (hi - lo);

            const auto candidate = trial_at(root);

            if (std::abs(candidate.residual) <=
                1.0e-8 * energy_scale) {
                roots.push_back(root);
            }
        }

        previous_pressure = pressure;
        previous = current;
    }

    if (roots.empty()) {
        throw std::runtime_error(
            "No admissible vaneless diffuser root found");
    }

    // For the reference model, select the highest-pressure
    // root. This is a branch-selection convention, not
    // a validated stability or choking criterion.
    const double p3 =
        *std::max_element(roots.begin(), roots.end());

    const auto outlet = trial_at(p3);

    const FlowStation station3{
        outlet.state,
        {outlet.meridional_velocity, whirl3}
    };

    const auto stagnation3 =
        evaluate_stagnation_state(model, station3);

    const double vm3 = outlet.meridional_velocity;

    return {
        state2,
        outlet.state,

        inlet.meridional_velocity,
        inlet.tangential_velocity,
        vm3,
        whirl3,

        area3,
        std::hypot(vm3, whirl3),

        p02,
        stagnation3.stagnation_state.p,

        outlet.state.rho * area3 * vm3 -
            inlet.mass_flow,

        stagnation3.stagnation_enthalpy - h02,

        outlet.state.s - s2,

        geometry.outlet_radius * whirl3 -
            geometry.inlet_radius *
                inlet.tangential_velocity,

        static_cast<int>(roots.size()),
        p3
    };
}

// Adapter connecting the existing impeller operating-point
// solution to the fixed-geometry diffuser.
//
// The existing impeller solver is not modified.

inline VanelessDiffuserResult
evaluate_vaneless_diffuser_from_operating_point(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& impeller_geometry,
    const OperatingPointRequest& request,
    const OperatingPointResult& impeller_solution,
    const VanelessDiffuserGeometry& diffuser_geometry)
{
    if (!(impeller_solution.converged &&
          impeller_solution.independently_validated &&
          impeller_solution.basic_physics_valid)) {
        throw std::invalid_argument(
            "Impeller solution is not independently verified");
    }

    const auto state2 = evaluate(
        model,
        impeller_solution.outlet_temperature,
        impeller_solution.outlet_pressure);

    const double area2 = effective_flow_area(
        impeller_exit_area(impeller_geometry.exit),
        impeller_geometry.exit.blockage);

    const double vm2 =
        request.operating.mass_flow /
        (state2.rho * area2);

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

    const double r2 =
        impeller_geometry.exit.diameter / 2.0;

    if (std::abs(diffuser_geometry.inlet_radius - r2) >
        1.0e-12 * std::max(1.0, r2)) {
        throw std::invalid_argument(
            "Diffuser inlet radius differs from impeller exit");
    }

    return evaluate_vaneless_diffuser(
        model,
        diffuser_geometry,
        {state2, vm2, whirl2, request.operating.mass_flow});
}

} // namespace turbocompgpu
