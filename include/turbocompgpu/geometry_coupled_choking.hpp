#pragma once

#include <turbocompgpu/inlet_coupled_choking.hpp>
#include <turbocompgpu/throat_geometry.hpp>
#include <turbocompgpu/throat_velocity.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

// Stage 6C: inlet-coupled, geometry-derived, blade-aligned
// rotating-throat choking.
//
// Assumptions:
//   * steady, adiabatic, isentropic flow
//   * blade-aligned relative velocity at the throat
//   * zero throat-local tangential relative velocity
//   * locally one-dimensional sonic condition Wn = a
//   * effective throat area supplied by Stage 6A
//
// The geometric throat is prescribed; this solver does not
// locate the minimum-area section of a 3D blade passage.

struct GeometryCoupledChokingInputs {
    double inlet_temperature;
    double inlet_pressure;
    double rotational_speed_rpm;
    double inlet_prewhirl;

    ThroatGeometry throat_geometry;
};

struct GeometryCoupledChokingResult {
    InletCoupledChokingResult coupled;
    ThroatGeometryResult geometry;

    ImpellerRelativeVelocity throat_relative_velocity;
    ThroatLocalVelocity throat_local_velocity;

    double throat_blade_speed;
    double throat_rothalpy;

    double rothalpy_residual;
    double sonic_residual;
    double mass_flow_residual;
    double continuity_residual;
};

inline GeometryCoupledChokingResult
evaluate_geometry_coupled_choking(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& impeller_geometry,
    const GeometryCoupledChokingInputs& input)
{
    if (!(std::isfinite(input.inlet_temperature) &&
          input.inlet_temperature > 0.0 &&
          std::isfinite(input.inlet_pressure) &&
          input.inlet_pressure > 0.0 &&
          std::isfinite(input.rotational_speed_rpm) &&
          input.rotational_speed_rpm >= 0.0 &&
          std::isfinite(input.inlet_prewhirl))) {
        throw std::invalid_argument(
            "Invalid geometry-coupled choking inputs");
    }

    const auto throat_geometry =
        evaluate_throat_geometry(input.throat_geometry);

    // In blade-aligned flow, the relative velocity tangent
    // to the local throat section is zero.
    //
    // This is NOT the impeller circumferential component.
    const InletCoupledChokingInputs coupled_input{
        input.inlet_temperature,
        input.inlet_pressure,
        input.rotational_speed_rpm,
        input.inlet_prewhirl,

        input.throat_geometry.radius,
        0.0,
        throat_geometry.effective_area
    };

    const auto coupled =
        evaluate_inlet_coupled_choking(
            model,
            derivatives,
            impeller_geometry,
            coupled_input);

    const auto basis =
        make_throat_velocity_basis(
            input.throat_geometry.blade_angle_deg);

    const double a =
        coupled.throat.critical.sound_speed;

    // Blade-aligned sonic throat:
    // Wn = a, Wtau = 0.
    const ThroatLocalVelocity local{
        a, 0.0
    };

    const auto velocity =
        from_throat_local(basis, local);

    const auto recovered =
        to_throat_local(basis, velocity);

    const double omega =
        2.0 * std::numbers::pi *
        input.rotational_speed_rpm / 60.0;

    const double Ut =
        omega * input.throat_geometry.radius;

    const double throat_rothalpy =
        coupled.throat.critical.enthalpy +
        0.5 * relative_speed_squared(velocity) -
        0.5 * Ut * Ut;

    const double rothalpy_residual =
        throat_rothalpy -
        coupled.inlet_rothalpy;

    const double sonic_residual =
        recovered.normal - a;

    const double mass_flow_residual =
        coupled.critical_mass_flow -
        throat_geometry.effective_area *
            coupled.throat.critical.density *
            recovered.normal;

    const double inlet_area =
        effective_flow_area(
            inlet_annulus_area(impeller_geometry.inlet),
            impeller_geometry.inlet.blockage);

    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    if (!inlet) {
        throw std::runtime_error(
            "Geometry-coupled inlet EOS evaluation failed");
    }

    const double inlet_mass_flow =
        inlet.value.rho *
        inlet_area *
        coupled.inlet_meridional_velocity;

    const double continuity_residual =
        inlet_mass_flow -
        coupled.critical_mass_flow;

    if (!(std::isfinite(throat_rothalpy) &&
          std::isfinite(rothalpy_residual) &&
          std::isfinite(sonic_residual) &&
          std::isfinite(mass_flow_residual) &&
          std::isfinite(continuity_residual))) {
        throw std::runtime_error(
            "Nonfinite geometry-coupled choking result");
    }

    return {
        coupled,
        throat_geometry,
        velocity,
        recovered,
        Ut,
        throat_rothalpy,
        rothalpy_residual,
        sonic_residual,
        mass_flow_residual,
        continuity_residual
    };
}

} // namespace turbocompgpu
