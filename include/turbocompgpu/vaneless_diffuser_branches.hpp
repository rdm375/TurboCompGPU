#pragma once

#include <turbocompgpu/vaneless_diffuser.hpp>
#include <turbocompgpu/isentropic_sound_speed.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace turbocompgpu {

// Stage 8B.4c: branch-aware lossless diffuser.
//
// The branch classification is based on the meridional
// Mach number, not static pressure recovery.
//
// The critical condition is Vm = a for the fixed-whirl,
// fixed-area outlet characteristic.
//
// This is an algebraic critical-capacity condition, not
// a complete multidimensional choking model.

enum class DiffuserBranch {
    Subsonic,
    Supersonic
};

enum class DiffuserBranchStatus {
    Converged,
    Critical,
    AboveCriticalCapacity,
    BranchUnavailable,
    ThermodynamicDomainLimit
};

struct DiffuserBranchPoint {
    thermomodel::State state;

    double pressure = 0.0;
    double meridional_velocity = 0.0;
    double tangential_velocity = 0.0;
    double sound_speed = 0.0;
    double meridional_mach = 0.0;
    double mass_capacity = 0.0;
};

struct DiffuserBranchResult {
    DiffuserBranchStatus status =
        DiffuserBranchStatus::BranchUnavailable;

    DiffuserBranch branch = DiffuserBranch::Subsonic;

    DiffuserBranchPoint outlet;
    DiffuserBranchPoint critical;

    double requested_mass_flow = 0.0;
    double critical_mass_flow = 0.0;

    double mass_residual = 0.0;
    double energy_residual = 0.0;
    double entropy_residual = 0.0;
    double angular_momentum_residual = 0.0;

    int sonic_iterations = 0;
    int pressure_iterations = 0;

    double final_pressure_bracket = 0.0;
    double relative_pressure_bracket = 0.0;

    bool low_pressure_domain_limited = false;
    double lowest_valid_pressure = 0.0;
};

inline DiffuserBranchResult
evaluate_vaneless_diffuser_branch(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const VanelessDiffuserGeometry& geometry,
    const VanelessDiffuserInlet& inlet,
    DiffuserBranch branch)
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
            "Invalid diffuser branch geometry or inlet");
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
        {
            inlet.meridional_velocity,
            inlet.tangential_velocity
        }
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, station2);

    const double h02 =
        stagnation2.stagnation_enthalpy;

    const double s2 = state2.s;
    const double p02 =
        stagnation2.stagnation_state.p;

    const auto limits = model.metadata();

    // Restrict the search to the valid EOS pressure range.
    const double p_min =
        std::max(
            limits.pmin * (1.0 + 1.0e-8),
            state2.p * 0.01);

    if (!(p_min < state2.p && p02 > state2.p)) {
        throw std::runtime_error(
            "Invalid diffuser characteristic pressure bounds");
    }

    auto point_at = [&](double pressure)
        -> DiffuserBranchPoint
    {
        const auto inversion = invert_temperature(
            model,
            pressure,
            s2,
            Property::Entropy,
            state2.T);

        const auto& state = inversion.state;

        const double vm_squared =
            2.0 * (h02 - state.h) -
            whirl3 * whirl3;

        if (!(vm_squared >= 0.0)) {
            throw std::domain_error(
                "No real meridional velocity");
        }

        const double vm = std::sqrt(vm_squared);

        const auto acoustic =
            evaluate_isentropic_sound_speed(
                derivatives, state.T, pressure);

        return {
            state,
            pressure,
            vm,
            whirl3,
            acoustic.sound_speed,
            vm / acoustic.sound_speed,
            state.rho * area3 * vm
        };
    };

    DiffuserBranchResult result;
    result.branch = branch;
    result.requested_mass_flow = inlet.mass_flow;

    // Establish a valid low-pressure bound by stepping
    // down from the known-valid impeller exit state.
    //
    // The final interval is not extrapolated across an
    // invalid thermodynamic state.

    double low_pressure = state2.p;
    DiffuserBranchPoint low = point_at(low_pressure);

    bool domain_limited = false;

    constexpr int scan_intervals = 512;

    for (int i = 1; i <= scan_intervals; ++i) {
        const double fraction =
            static_cast<double>(i) / scan_intervals;

        const double candidate_pressure =
            state2.p *
            std::exp(
                fraction * std::log(p_min / state2.p));

        try {
            const auto candidate =
                point_at(candidate_pressure);

            low_pressure = candidate_pressure;
            low = candidate;
        }
        catch (const std::exception&) {
            domain_limited = true;
            break;
        }
    }

    result.low_pressure_domain_limited = domain_limited;
    result.lowest_valid_pressure = low_pressure;

    // Locate the sonic point using the Mach residual.
    //
    // At the lower bound we require supersonic flow.
    // At the impeller exit, or at a higher pressure,
    // we require subsonic flow.

    double sonic_lo = low_pressure;
    double sonic_hi = state2.p;

    auto sonic_low = low;
    auto sonic_high = point_at(sonic_hi);

    if (sonic_high.meridional_mach >= 1.0) {
        // Search upward in pressure for a subsonic state.
        constexpr int upper_intervals = 512;

        bool found = false;

        for (int i = 1; i <= upper_intervals; ++i) {
            const double fraction =
                static_cast<double>(i) / upper_intervals;

            const double pressure =
                state2.p + fraction * (p02 - state2.p);

            try {
                const auto candidate = point_at(pressure);

                if (candidate.meridional_mach < 1.0) {
                    sonic_hi = pressure;
                    sonic_high = candidate;
                    found = true;
                    break;
                }
            }
            catch (const std::exception&) {
                break;
            }
        }

        if (!found)
            throw std::runtime_error(
                "Could not bracket subsonic side");
    }

    if (sonic_low.meridional_mach <= 1.0) {
        throw std::runtime_error(
            "Could not bracket supersonic side");
    }

    for (int iteration = 0; iteration < 120; ++iteration) {
        const double mid =
            sonic_lo + 0.5 * (sonic_hi - sonic_lo);

        if (mid == sonic_lo || mid == sonic_hi)
            break;

        const auto point = point_at(mid);

        ++result.sonic_iterations;

        if (point.meridional_mach > 1.0) {
            sonic_lo = mid;
        } else {
            sonic_hi = mid;
        }
    }

    result.critical =
        point_at(sonic_lo + 0.5 * (sonic_hi - sonic_lo));

    result.critical_mass_flow =
        result.critical.mass_capacity;

    const double mass_tolerance =
        1.0e-10 *
        std::max(1.0, result.critical_mass_flow);

    if (inlet.mass_flow >
        result.critical_mass_flow + mass_tolerance) {
        result.status =
            DiffuserBranchStatus::AboveCriticalCapacity;
        return result;
    }

    if (std::abs(
            inlet.mass_flow - result.critical_mass_flow) <=
        mass_tolerance) {
        result.status = DiffuserBranchStatus::Critical;
        result.outlet = result.critical;
    } else {
        // For the subsonic branch, capacity decreases as
        // pressure rises above the sonic pressure.
        //
        // For the supersonic branch, capacity decreases
        // as pressure falls below the sonic pressure.

        double lo = 0.0;
        double hi = 0.0;

        if (branch == DiffuserBranch::Supersonic) {
            lo = low_pressure;
            hi = result.critical.pressure;

            if (low.mass_capacity > inlet.mass_flow) {
                result.status = domain_limited
                    ? DiffuserBranchStatus::ThermodynamicDomainLimit
                    : DiffuserBranchStatus::BranchUnavailable;
                return result;
            }
        } else {
            lo = result.critical.pressure;

            // Locate the zero-meridional-velocity limit
            // by scanning upward toward stagnation pressure.
            //
            // The final admissible state must have capacity
            // below the requested mass flow.

            hi = state2.p;

            bool found = false;

            for (int i = 0; i <= scan_intervals; ++i) {
                const double fraction =
                    static_cast<double>(i) / scan_intervals;

                const double pressure =
                    state2.p + fraction * (p02 - state2.p);

                try {
                    const auto candidate =
                        point_at(pressure);

                    if (candidate.mass_capacity <
                        inlet.mass_flow) {
                        hi = pressure;
                        found = true;
                        break;
                    }
                }
                catch (const std::exception&) {
                    break;
                }
            }

            if (!found) {
                result.status =
                    DiffuserBranchStatus::BranchUnavailable;
                return result;
            }
        }

        // Both branches are now bracketed in pressure.
        // The capacity residual changes sign monotonically
        // on the selected branch.

        double flo =
            point_at(lo).mass_capacity - inlet.mass_flow;

        double fhi =
            point_at(hi).mass_capacity - inlet.mass_flow;

        if (flo * fhi > 0.0) {
            result.status =
                DiffuserBranchStatus::BranchUnavailable;
            return result;
        }

        for (int iteration = 0; iteration < 120; ++iteration) {
            const double mid =
                lo + 0.5 * (hi - lo);

            if (mid == lo || mid == hi)
                break;

            const auto point = point_at(mid);

            ++result.pressure_iterations;

            const double fm =
                point.mass_capacity - inlet.mass_flow;

            if ((fm < 0.0) == (flo < 0.0)) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
                fhi = fm;
            }
        }

        result.final_pressure_bracket = hi - lo;

        result.relative_pressure_bracket =
            result.final_pressure_bracket /
            std::max(
                1.0,
                std::abs(lo + 0.5 * (hi - lo)));

        result.outlet =
            point_at(lo + 0.5 * (hi - lo));

        result.status = DiffuserBranchStatus::Converged;
    }

    // Independent conservation diagnostics.

    const auto& outlet = result.outlet;

    result.mass_residual =
        outlet.mass_capacity - inlet.mass_flow;

    result.energy_residual =
        outlet.state.h +
        0.5 * (
            outlet.meridional_velocity *
                outlet.meridional_velocity +
            outlet.tangential_velocity *
                outlet.tangential_velocity) -
        h02;

    result.entropy_residual =
        outlet.state.s - s2;

    result.angular_momentum_residual =
        geometry.outlet_radius *
            outlet.tangential_velocity -
        geometry.inlet_radius *
            inlet.tangential_velocity;

    return result;
}

} // namespace turbocompgpu
