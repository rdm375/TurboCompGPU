#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace turbocompgpu {

// Thermodynamic data needed to integrate a constant-polytropic-efficiency
// compression path in pressure coordinates.
struct PolytropicPathProperties {
    double density_kg_m3;
    double dh_dp_T_j_kg_pa;
    double dh_dT_p_j_kg_k;
};

struct PolytropicEfficiencyResult {
    double polytropic_head_j_kg;
    double polytropic_efficiency;
    double endpoint_temperature_residual_k;
    std::size_t pressure_steps;
    std::size_t efficiency_iterations;
};

namespace polytropic_efficiency_detail {

struct PathState {
    double temperature;
    double ideal_head;
};

template<class EvaluateProperties>
PathState integrate_constant_efficiency_path(
    double inlet_temperature,
    double inlet_pressure,
    double outlet_pressure,
    double efficiency,
    std::size_t steps,
    EvaluateProperties&& evaluate_properties)
{
    const double log_pressure_ratio =
        std::log(outlet_pressure / inlet_pressure);
    const double dx = log_pressure_ratio / static_cast<double>(steps);

    auto rhs = [&](double x, const PathState& y) {
        const double pressure = inlet_pressure * std::exp(x);
        const auto properties = evaluate_properties(y.temperature, pressure);
        if (!(std::isfinite(properties.density_kg_m3) &&
              properties.density_kg_m3 > 0.0 &&
              std::isfinite(properties.dh_dp_T_j_kg_pa) &&
              std::isfinite(properties.dh_dT_p_j_kg_k) &&
              properties.dh_dT_p_j_kg_k > 0.0)) {
            throw std::runtime_error(
                "Invalid EOS properties in polytropic path integration");
        }

        const double specific_volume = 1.0 / properties.density_kg_m3;
        // dh = T ds + v dp and eta_p = v dp / dh, so for constant eta_p:
        // dh/dp = v/eta_p. Expand dh/dp along T(p) with EOS derivatives.
        const double dT_dp =
            (specific_volume / efficiency -
             properties.dh_dp_T_j_kg_pa) /
            properties.dh_dT_p_j_kg_k;
        return PathState{
            pressure * dT_dp,
            pressure * specific_volume
        };
    };

    PathState y{inlet_temperature, 0.0};
    double x = 0.0;
    for (std::size_t i = 0; i < steps; ++i) {
        const auto k1 = rhs(x, y);
        const auto k2 = rhs(x + 0.5 * dx,
            {y.temperature + 0.5 * dx * k1.temperature,
             y.ideal_head + 0.5 * dx * k1.ideal_head});
        const auto k3 = rhs(x + 0.5 * dx,
            {y.temperature + 0.5 * dx * k2.temperature,
             y.ideal_head + 0.5 * dx * k2.ideal_head});
        const auto k4 = rhs(x + dx,
            {y.temperature + dx * k3.temperature,
             y.ideal_head + dx * k3.ideal_head});
        y.temperature += dx / 6.0 *
            (k1.temperature + 2.0 * k2.temperature +
             2.0 * k3.temperature + k4.temperature);
        y.ideal_head += dx / 6.0 *
            (k1.ideal_head + 2.0 * k2.ideal_head +
             2.0 * k3.ideal_head + k4.ideal_head);
        x += dx;
        if (!(std::isfinite(y.temperature) && y.temperature > 0.0 &&
              std::isfinite(y.ideal_head) && y.ideal_head > 0.0)) {
            throw std::runtime_error(
                "Polytropic pressure path left the physical domain");
        }
    }
    return y;
}

} // namespace polytropic_efficiency_detail

// Infer eta_p by shooting a constant-efficiency real-gas path from the inlet
// until its outlet temperature matches the supplied actual outlet state.
// The caller supplies EOS derivatives on a mass-specific basis. This is an
// independent diagnostic; it does not change a stage solve or eta_tt.
template<class EvaluateProperties>
PolytropicEfficiencyResult calculate_polytropic_efficiency(
    double inlet_temperature_k,
    double inlet_pressure_pa,
    double inlet_enthalpy_j_kg,
    double outlet_temperature_k,
    double outlet_pressure_pa,
    double outlet_enthalpy_j_kg,
    EvaluateProperties&& evaluate_properties,
    std::size_t pressure_steps = 128,
    std::size_t efficiency_iterations = 48)
{
    if (!(std::isfinite(inlet_temperature_k) && inlet_temperature_k > 0.0 &&
          std::isfinite(outlet_temperature_k) && outlet_temperature_k > 0.0 &&
          std::isfinite(inlet_pressure_pa) && inlet_pressure_pa > 0.0 &&
          std::isfinite(outlet_pressure_pa) &&
          outlet_pressure_pa > inlet_pressure_pa &&
          std::isfinite(inlet_enthalpy_j_kg) &&
          std::isfinite(outlet_enthalpy_j_kg) &&
          outlet_enthalpy_j_kg > inlet_enthalpy_j_kg &&
          pressure_steps >= 4 && efficiency_iterations > 0)) {
        throw std::invalid_argument("Invalid polytropic efficiency inputs");
    }

    auto endpoint_residual = [&](double efficiency) {
        return polytropic_efficiency_detail::integrate_constant_efficiency_path(
            inlet_temperature_k, inlet_pressure_pa, outlet_pressure_pa,
            efficiency, pressure_steps, evaluate_properties).temperature -
            outlet_temperature_k;
    };

    double upper = 1.0;
    const double upper_residual = endpoint_residual(upper);
    if (!std::isfinite(upper_residual))
        throw std::runtime_error("Non-finite polytropic path residual");
    if (upper_residual > 1e-8) {
        throw std::runtime_error(
            "Outlet temperature is below the EOS isentropic limit");
    }
    double lower = upper;
    bool bracketed = false;
    while (lower > 0.02) {
        const double candidate = std::max(0.02, 0.5 * lower);
        double candidate_residual = 0.0;
        try {
            candidate_residual = endpoint_residual(candidate);
        } catch (const std::exception&) {
            // Very low eta can leave the EOS validity envelope before the
            // physical root. Treat that as the high-temperature side of the
            // bracket; bisection will move back into the model domain.
            candidate_residual = INFINITY;
        }
        if (candidate_residual >= 0.0) {
            lower = candidate;
            bracketed = true;
            break;
        }
        upper = candidate;
        lower = candidate;
    }
    if (!bracketed) {
        throw std::runtime_error(
            "Outlet temperature implies polytropic efficiency below 0.02");
    }

    // `upper` is the higher eta (lower-temperature) side and `lower` is the
    // lower eta (higher-temperature) side of the endpoint root.
    double eta_hi = upper;
    double eta_lo = lower;
    for (std::size_t i = 0; i < efficiency_iterations; ++i) {
        const double middle = 0.5 * (eta_lo + eta_hi);
        double residual = 0.0;
        try {
            residual = endpoint_residual(middle);
        } catch (const std::exception&) {
            residual = INFINITY;
        }
        if (residual >= 0.0)
            eta_lo = middle;
        else
            eta_hi = middle;
    }

    const double eta = 0.5 * (eta_lo + eta_hi);
    const auto path =
        polytropic_efficiency_detail::integrate_constant_efficiency_path(
            inlet_temperature_k, inlet_pressure_pa, outlet_pressure_pa,
            eta, pressure_steps, evaluate_properties);
    const double actual_enthalpy_rise =
        outlet_enthalpy_j_kg - inlet_enthalpy_j_kg;
    const double efficiency_from_head = path.ideal_head / actual_enthalpy_rise;
    if (!(std::isfinite(efficiency_from_head) &&
          efficiency_from_head > 0.0 && efficiency_from_head <= 1.0 + 1e-8)) {
        throw std::runtime_error(
            "Integrated polytropic head gives an invalid efficiency");
    }

    return {path.ideal_head, efficiency_from_head,
            path.temperature - outlet_temperature_k,
            pressure_steps, efficiency_iterations};
}

} // namespace turbocompgpu
