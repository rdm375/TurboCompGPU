#pragma once

#include <turbocompgpu/map_coordinates.hpp>
#include <turbocompgpu/operating_point.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct CompressorMapPoint {
    double mass_flow;
    double rotational_speed_rpm;

    double corrected_mass_flow;
    double corrected_speed_rpm;

    double inlet_total_temperature;
    double inlet_total_pressure;
    double outlet_total_temperature;
    double outlet_total_pressure;

    double static_pressure_ratio;
    double total_pressure_ratio;

    double specific_work;
    double shaft_power;
    double torque;
    double entropy_generation;

    double full_residual_norm;
    double inlet_stagnation_enthalpy_residual;
    double inlet_stagnation_entropy_residual;
    double outlet_stagnation_enthalpy_residual;
    double outlet_stagnation_entropy_residual;

    OperatingLimitStatus operating_limits;
};

// Postprocessing only: the nonlinear solver is not called.
//
// Requires an already converged, independently validated,
// basic-physics-valid operating point.
//
// Uses the existing full-system reconstruction and its
// EOS-consistent inlet/outlet stagnation states.
inline CompressorMapPoint evaluate_compressor_map_point(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointResult& point,
    const NonlinearResidualScales& scales,
    const MapReferenceConditions& reference,
    double consistency_tolerance = 1.0e-8)
{
    if (!point.converged ||
        !point.independently_validated ||
        !point.basic_physics_valid) {
        throw std::invalid_argument(
            "Compressor map requires a validated operating point");
    }

    if (!(std::isfinite(consistency_tolerance) &&
          consistency_tolerance > 0.0)) {
        throw std::invalid_argument(
            "Invalid compressor-map consistency tolerance");
    }

    const NonlinearUnknowns unknowns {
        point.outlet_temperature,
        std::log(point.outlet_pressure / scales.pressure)
    };

    const auto full = evaluate_nonlinear_residual(
        model,
        geometry,
        point.operating,
        unknowns,
        scales
    );

    if (!(std::isfinite(full.residual_norm) &&
          full.residual_norm <= consistency_tolerance)) {
        throw std::runtime_error(
            "Compressor-map full-system validation failed");
    }

    const auto& system = full.full_system.system;

    const auto& inlet_total =
        system.inlet_stagnation.stagnation_state;

    const auto& outlet_total =
        system.outlet_stagnation.stagnation_state;

    const auto close = [consistency_tolerance](
        double a, double b) {
        return std::isfinite(a) &&
               std::isfinite(b) &&
               std::abs(a - b) <=
                   consistency_tolerance *
                   std::max({1.0, std::abs(a), std::abs(b)});
    };

    if (!close(full.outlet_pressure, point.outlet_pressure) ||
        !close(full.euler_specific_work, point.specific_work) ||
        !close(full.shaft_torque, point.torque)) {
        throw std::runtime_error(
            "Compressor-map reconstruction disagrees with operating point");
    }

    const auto corrected = corrected_map_coordinates(
        point.operating.mass_flow,
        point.operating.rotational_speed_rpm,
        {inlet_total.T, inlet_total.p},
        reference
    );

    const auto static_ratio = map_pressure_ratio(
        system.inlet_static.p,
        system.outlet_static.p,
        MapPressureConvention::StaticToStatic
    );

    const auto total_ratio = map_pressure_ratio(
        inlet_total.p,
        outlet_total.p,
        MapPressureConvention::TotalToTotal
    );

    return {
        point.operating.mass_flow,
        point.operating.rotational_speed_rpm,

        corrected.corrected_mass_flow,
        corrected.corrected_speed_rpm,

        inlet_total.T,
        inlet_total.p,
        outlet_total.T,
        outlet_total.p,

        static_ratio.value,
        total_ratio.value,

        point.specific_work,
        point.shaft_power,
        point.torque,
        point.entropy_generation,

        full.residual_norm,

        system.inlet_stagnation.enthalpy_residual,
        system.inlet_stagnation.entropy_residual,
        system.outlet_stagnation.enthalpy_residual,
        system.outlet_stagnation.entropy_residual,

        point.operating_limits
    };
}

} // namespace turbocompgpu
