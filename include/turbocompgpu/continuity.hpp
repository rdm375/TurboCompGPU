#pragma once

#include <turbocompgpu/impeller.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct InletAnnulus {
    double hub_diameter;       // m
    double shroud_diameter;    // m
    double blockage;           // fraction [0,1)
};

struct ImpellerExitPassage {
    double diameter;           // m
    double width;              // m
    double blockage;           // fraction [0,1)
};

struct ContinuityOperatingPoint {
    double mass_flow;          // kg/s
    double inlet_density;      // kg/m^3
    double exit_density;       // kg/m^3
};

struct ContinuityResult {
    double inlet_geometric_area;  // m^2
    double inlet_effective_area;  // m^2
    double exit_geometric_area;   // m^2
    double exit_effective_area;   // m^2
    double inlet_meridional_velocity; // m/s
    double exit_meridional_velocity;  // m/s
    double inlet_mass_flow;      // kg/s
    double exit_mass_flow;       // kg/s
};

inline double inlet_annulus_area(
    const InletAnnulus& inlet)
{
    const double Dh = inlet.hub_diameter;
    const double Ds = inlet.shroud_diameter;

    if (!std::isfinite(Dh) ||
        !std::isfinite(Ds) ||
        Dh < 0.0 ||
        Ds <= Dh) {
        throw std::invalid_argument(
            "Invalid inlet hub/shroud diameters");
    }

    return std::numbers::pi / 4.0 *
           (Ds - Dh) * (Ds + Dh);
}

inline double impeller_exit_area(
    const ImpellerExitPassage& exit)
{
    if (!std::isfinite(exit.diameter) ||
        exit.diameter <= 0.0 ||
        !std::isfinite(exit.width) ||
        exit.width <= 0.0) {
        throw std::invalid_argument(
            "Invalid impeller exit diameter or width");
    }

    return std::numbers::pi *
           exit.diameter * exit.width;
}

inline double effective_flow_area(
    double geometric_area,
    double blockage)
{
    if (!std::isfinite(geometric_area) ||
        geometric_area <= 0.0 ||
        !std::isfinite(blockage) ||
        blockage < 0.0 ||
        blockage >= 1.0) {
        throw std::invalid_argument(
            "Invalid flow area or blockage fraction");
    }

    const double area =
        geometric_area * (1.0 - blockage);

    if (!std::isfinite(area) || area <= 0.0) {
        throw std::runtime_error(
            "Nonpositive effective flow area");
    }

    return area;
}

inline double meridional_velocity(
    double mass_flow,
    double density,
    double effective_area)
{
    if (!std::isfinite(mass_flow) ||
        mass_flow < 0.0 ||
        !std::isfinite(density) ||
        density <= 0.0 ||
        !std::isfinite(effective_area) ||
        effective_area <= 0.0) {
        throw std::invalid_argument(
            "Invalid mass flow, density, or area");
    }

    const double velocity =
        (mass_flow / density) / effective_area;

    if (!std::isfinite(velocity)) {
        throw std::runtime_error(
            "Non-finite meridional velocity");
    }

    return velocity;
}

inline ContinuityResult evaluate_continuity(
    const InletAnnulus& inlet,
    const ImpellerExitPassage& exit,
    const ContinuityOperatingPoint& operating_point)
{
    const double A1g = inlet_annulus_area(inlet);
    const double A2g = impeller_exit_area(exit);

    const double A1 =
        effective_flow_area(A1g, inlet.blockage);

    const double A2 =
        effective_flow_area(A2g, exit.blockage);

    const double Vm1 =
        meridional_velocity(
            operating_point.mass_flow,
            operating_point.inlet_density,
            A1);

    const double Vm2 =
        meridional_velocity(
            operating_point.mass_flow,
            operating_point.exit_density,
            A2);

    const double m1 =
        (operating_point.inlet_density * A1) * Vm1;

    const double m2 =
        (operating_point.exit_density * A2) * Vm2;

    if (!std::isfinite(m1) || !std::isfinite(m2)) {
        throw std::runtime_error(
            "Non-finite reconstructed mass flow");
    }

    return {
        A1g, A1,
        A2g, A2,
        Vm1, Vm2,
        m1, m2
    };
}

// Adapter: use continuity-derived velocities in the
// previously validated aerodynamic model.
inline ImpellerResult evaluate_impeller_from_continuity(
    const ImpellerGeometry& geometry,
    const ImpellerOperatingPoint& operating_point,
    const ContinuityResult& continuity)
{
    ImpellerOperatingPoint derived = operating_point;

    derived.inlet_meridional_velocity =
        continuity.inlet_meridional_velocity;

    derived.exit_meridional_velocity =
        continuity.exit_meridional_velocity;

    return evaluate_impeller(geometry, derived);
}

} // namespace turbocompgpu
