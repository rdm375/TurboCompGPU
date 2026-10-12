#pragma once

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// User-specified reference conditions.
// No particular industry reference convention is assumed.
struct MapReferenceConditions {
    double temperature;
    double pressure;
};

// Physical operating coordinates.
//
// Temperature and pressure must be inlet TOTAL conditions
// when conventional corrected coordinates are requested.
struct MapInletConditions {
    double total_temperature;
    double total_pressure;
};

struct CorrectedMapCoordinates {
    double corrected_mass_flow;
    double corrected_speed_rpm;

    double temperature_ratio;
    double pressure_ratio;
};

inline CorrectedMapCoordinates corrected_map_coordinates(
    double mass_flow,
    double rotational_speed_rpm,
    const MapInletConditions& inlet,
    const MapReferenceConditions& reference)
{
    const auto positive_finite = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };

    if (!positive_finite(mass_flow) ||
        !positive_finite(rotational_speed_rpm) ||
        !positive_finite(inlet.total_temperature) ||
        !positive_finite(inlet.total_pressure) ||
        !positive_finite(reference.temperature) ||
        !positive_finite(reference.pressure)) {
        throw std::invalid_argument(
            "Map coordinates require positive finite inputs");
    }

    const double temperature_ratio =
        inlet.total_temperature / reference.temperature;

    const double pressure_ratio =
        inlet.total_pressure / reference.pressure;

    const double corrected_mass_flow =
        mass_flow * std::sqrt(temperature_ratio) /
        pressure_ratio;

    const double corrected_speed_rpm =
        rotational_speed_rpm / std::sqrt(temperature_ratio);

    if (!positive_finite(corrected_mass_flow) ||
        !positive_finite(corrected_speed_rpm)) {
        throw std::overflow_error(
            "Corrected map coordinates are not finite");
    }

    return {
        corrected_mass_flow,
        corrected_speed_rpm,
        temperature_ratio,
        pressure_ratio
    };
}

// Explicit pressure-ratio conventions.
//
// M1.5 reports thermodynamic station-pressure ratios.
// Do not silently reinterpret those as total-to-total ratios.
enum class MapPressureConvention {
    StaticToStatic,
    TotalToTotal
};

struct MapPressureRatio {
    double value;
    MapPressureConvention convention;
};

inline MapPressureRatio map_pressure_ratio(
    double inlet_pressure,
    double outlet_pressure,
    MapPressureConvention convention)
{
    if (!(std::isfinite(inlet_pressure) &&
          inlet_pressure > 0.0 &&
          std::isfinite(outlet_pressure) &&
          outlet_pressure > 0.0)) {
        throw std::invalid_argument(
            "Map pressures must be positive and finite");
    }

    const double ratio = outlet_pressure / inlet_pressure;

    if (!std::isfinite(ratio) || ratio <= 0.0) {
        throw std::overflow_error("Invalid map pressure ratio");
    }

    return {ratio, convention};
}

} // namespace turbocompgpu
