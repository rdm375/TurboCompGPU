#pragma once

#include <turbocompgpu/operating_incidence.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

inline double zero_incidence_mass_flow(
    const ImpellerSystemGeometry& geometry,
    double inlet_density,
    double rotational_speed_rpm,
    double inlet_prewhirl,
    double blade_inlet_angle_deg)
{
    if (!std::isfinite(inlet_density) ||
        inlet_density <= 0.0 ||
        !std::isfinite(rotational_speed_rpm) ||
        rotational_speed_rpm <= 0.0 ||
        !std::isfinite(inlet_prewhirl) ||
        !std::isfinite(blade_inlet_angle_deg)) {
        throw std::invalid_argument(
            "Invalid zero-incidence locus inputs"
        );
    }

    const double area = effective_flow_area(
        inlet_annulus_area(geometry.inlet),
        geometry.inlet.blockage
    );

    const double D1 = geometry.impeller.inlet_diameter;
    const double D2 = geometry.impeller.exit_diameter;

    if (!std::isfinite(D1) ||
        !std::isfinite(D2) ||
        D1 <= 0.0 ||
        D1 >= D2) {
        throw std::invalid_argument(
            "Invalid impeller inlet diameter"
        );
    }

    const double U1 =
        std::numbers::pi * D1 *
        rotational_speed_rpm / 60.0;

    const double beta =
        blade_inlet_angle_deg *
        std::numbers::pi / 180.0;

    // This correlation uses the conventional
    // positive-meridional-flow inlet triangle.
    if (!(U1 > inlet_prewhirl) ||
        !(beta > 0.0) ||
        !(beta < 0.5 * std::numbers::pi)) {
        throw std::invalid_argument(
            "Zero-incidence locus outside supported regime"
        );
    }

    const double mass_flow =
        inlet_density * area *
        (U1 - inlet_prewhirl) *
        std::tan(beta);

    if (!std::isfinite(mass_flow) ||
        mass_flow <= 0.0) {
        throw std::runtime_error(
            "Invalid zero-incidence mass flow"
        );
    }

    return mass_flow;
}

} // namespace turbocompgpu
