#pragma once

#include <turbocompgpu/impeller_system.hpp>
#include <turbocompgpu/throat_velocity.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

// Stage 6E: inlet velocity triangle and prescribed
// inlet-to-throat kinematic transport.
//
// The initial transport closure preserves both relative
// velocity components. It is a baseline approximation,
// not a blade-passage momentum or turning model.

struct ThroatIncidenceInputs {
    double inlet_density;
    double mass_flow;
    double rotational_speed_rpm;
    double inlet_prewhirl;
    double throat_blade_angle_deg;
};

struct ThroatIncidenceResult {
    double inlet_area;
    double inlet_meridional_velocity;
    double inlet_blade_speed;

    ImpellerRelativeVelocity inlet_relative_velocity;
    ImpellerRelativeVelocity transported_relative_velocity;

    ThroatLocalVelocity transported_local_velocity;

    double incidence_angle_rad;
    double incidence_angle_deg;
};

inline ThroatIncidenceResult evaluate_throat_incidence(
    const ImpellerSystemGeometry& geometry,
    const ThroatIncidenceInputs& input)
{
    if (!(std::isfinite(input.inlet_density) &&
          input.inlet_density > 0.0 &&
          std::isfinite(input.mass_flow) &&
          input.mass_flow >= 0.0 &&
          std::isfinite(input.rotational_speed_rpm) &&
          input.rotational_speed_rpm >= 0.0 &&
          std::isfinite(input.inlet_prewhirl))) {
        throw std::invalid_argument(
            "Invalid throat-incidence inputs");
    }

    const double inlet_area =
        effective_flow_area(
            inlet_annulus_area(geometry.inlet),
            geometry.inlet.blockage);

    if (!(std::isfinite(inlet_area) && inlet_area > 0.0))
        throw std::invalid_argument("Invalid inlet flow area");

    const double omega =
        2.0 * std::numbers::pi *
        input.rotational_speed_rpm / 60.0;

    const double inlet_radius =
        geometry.impeller.inlet_diameter * 0.5;

    if (!(std::isfinite(inlet_radius) &&
          inlet_radius >= 0.0)) {
        throw std::invalid_argument("Invalid impeller inlet radius");
    }

    const double U1 = omega * inlet_radius;

    const double Vm1 =
        input.mass_flow / (input.inlet_density * inlet_area);

    const ImpellerRelativeVelocity inlet_relative{
        Vm1,
        input.inlet_prewhirl - U1
    };

    // Explicit baseline transport closure.
    const ImpellerRelativeVelocity transported =
        inlet_relative;

    const auto basis =
        make_throat_velocity_basis(
            input.throat_blade_angle_deg);

    const auto local =
        to_throat_local(basis, transported);

    // Signed angle relative to the throat-normal direction.
    const double angle =
        std::atan2(local.tangential, local.normal);

    if (!(std::isfinite(Vm1) &&
          std::isfinite(U1) &&
          std::isfinite(angle))) {
        throw std::runtime_error(
            "Nonfinite throat-incidence result");
    }

    return {
        inlet_area,
        Vm1,
        U1,
        inlet_relative,
        transported,
        local,
        angle,
        angle * 180.0 / std::numbers::pi
    };
}

} // namespace turbocompgpu
