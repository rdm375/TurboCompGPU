#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

// Idealized local blade-passage throat geometry.
//
// Conventions:
//   radius:              throat meanline radius [m]
//   passage_height:      local blade-span height [m]
//   blade_angle_deg:     blade angle from tangential direction [deg]
//   blade_thickness:     thickness normal to blade surface [m]
//   blockage:            additional fractional flow-area blockage
//
// Blade angle is restricted to (0,90] degrees.
//
// This is a local geometric approximation. It does not
// locate the minimum-area throat of a 3D impeller.

struct ThroatGeometry {
    int blade_count;
    double radius;
    double passage_height;
    double blade_angle_deg;
    double blade_thickness;
    double blockage;
};

struct ThroatGeometryResult {
    double circumferential_pitch;
    double projected_pitch;
    double open_width;

    double area_per_passage;
    double geometric_area;
    double effective_area;

    double blade_angle_rad;
    double open_area_fraction;
};

inline ThroatGeometryResult evaluate_throat_geometry(
    const ThroatGeometry& geometry)
{
    const int N = geometry.blade_count;
    const double r = geometry.radius;
    const double b = geometry.passage_height;
    const double beta_deg = geometry.blade_angle_deg;
    const double t = geometry.blade_thickness;
    const double B = geometry.blockage;

    if (N < 2) {
        throw std::invalid_argument(
            "Throat blade count must be at least two");
    }

    if (!(std::isfinite(r) && r > 0.0 &&
          std::isfinite(b) && b > 0.0 &&
          std::isfinite(beta_deg) &&
          beta_deg > 0.0 &&
          beta_deg <= 90.0 &&
          std::isfinite(t) && t >= 0.0 &&
          std::isfinite(B) &&
          B >= 0.0 && B < 1.0)) {
        throw std::invalid_argument(
            "Invalid throat geometry parameters");
    }

    const double beta =
        beta_deg * std::numbers::pi / 180.0;

    const double pitch =
        2.0 * std::numbers::pi * r /
        static_cast<double>(N);

    // Pitch projected onto the blade-normal direction.
    const double projected_pitch =
        pitch * std::sin(beta);

    const double open_width =
        projected_pitch - t;

    if (!(std::isfinite(open_width) &&
          open_width > 0.0)) {
        throw std::invalid_argument(
            "Blade thickness closes the throat passage");
    }

    const double area_per_passage =
        b * open_width;

    const double geometric_area =
        static_cast<double>(N) *
        area_per_passage;

    const double effective_area =
        geometric_area * (1.0 - B);

    const double open_fraction =
        open_width / projected_pitch;

    if (!(std::isfinite(effective_area) &&
          effective_area > 0.0 &&
          std::isfinite(open_fraction) &&
          open_fraction > 0.0 &&
          open_fraction <= 1.0)) {
        throw std::runtime_error(
            "Invalid calculated throat area");
    }

    return {
        pitch,
        projected_pitch,
        open_width,
        area_per_passage,
        geometric_area,
        effective_area,
        beta,
        open_fraction
    };
}

} // namespace turbocompgpu
