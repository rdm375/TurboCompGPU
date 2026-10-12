#pragma once

#include <turbocompgpu/impeller_system.hpp>
#include <turbocompgpu/incidence_loss.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct OperatingIncidenceConditions {
    double mass_flow;             // kg/s
    double rotational_speed_rpm;
    double inlet_density;         // kg/m^3
    double inlet_prewhirl;        // m/s
};

struct OperatingIncidenceResult {
    double inlet_effective_area;
    double inlet_meridional_velocity;
    double inlet_blade_speed;
    double inlet_relative_velocity;
    double inlet_relative_angle_deg;
    IncidenceLossResult loss;
};

inline OperatingIncidenceResult evaluate_operating_incidence_loss(
    const ImpellerSystemGeometry& geometry,
    const OperatingIncidenceConditions& operating,
    const IncidenceLossParameters& parameters)
{
    const double area = effective_flow_area(
        inlet_annulus_area(geometry.inlet),
        geometry.inlet.blockage
    );

    const double Vm1 = meridional_velocity(
        operating.mass_flow,
        operating.inlet_density,
        area
    );

    const double D1 = geometry.impeller.inlet_diameter;
    const double D2 = geometry.impeller.exit_diameter;

    if (!std::isfinite(D1) ||
        !std::isfinite(D2) ||
        D1 <= 0.0 ||
        D1 >= D2 ||
        !std::isfinite(operating.rotational_speed_rpm) ||
        operating.rotational_speed_rpm <= 0.0 ||
        !std::isfinite(operating.inlet_prewhirl)) {
        throw std::invalid_argument(
            "Invalid operating-incidence geometry or conditions"
        );
    }

    const double U1 =
        std::numbers::pi *
        D1 *
        operating.rotational_speed_rpm / 60.0;

    const double Wtheta1 =
        operating.inlet_prewhirl - U1;

    const double W1 = std::hypot(
        Vm1, Wtheta1
    );

    const double beta1 =
        std::atan2(Vm1, -Wtheta1) *
        180.0 / std::numbers::pi;

    const auto loss = evaluate_incidence_loss(
        beta1, W1, parameters
    );

    return {
        area,
        Vm1,
        U1,
        W1,
        beta1,
        loss
    };
}

inline double calibrate_inlet_blade_angle_deg(
    const ImpellerSystemGeometry& geometry,
    const OperatingIncidenceConditions& design)
{
    const auto result = evaluate_operating_incidence_loss(
        geometry,
        design,
        {0.0, 0.0}
    );

    return result.inlet_relative_angle_deg;
}

} // namespace turbocompgpu
