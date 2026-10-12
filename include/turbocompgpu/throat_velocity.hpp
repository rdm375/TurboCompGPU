#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace turbocompgpu {

struct ImpellerRelativeVelocity {
    double meridional;
    double circumferential;
};

struct ThroatLocalVelocity {
    double normal;
    double tangential;
};

struct ThroatVelocityBasis {
    double blade_angle_rad;

    // Positive throat-normal unit vector.
    double normal_meridional;
    double normal_circumferential;

    // Positive throat-tangential unit vector.
    double tangent_meridional;
    double tangent_circumferential;
};

inline ThroatVelocityBasis make_throat_velocity_basis(
    double blade_angle_deg)
{
    if (!std::isfinite(blade_angle_deg) ||
        blade_angle_deg <= 0.0 ||
        blade_angle_deg > 90.0) {
        throw std::invalid_argument(
            "Throat blade angle must be in (0,90] degrees");
    }

    const double beta =
        blade_angle_deg * std::numbers::pi / 180.0;

    const double s = std::sin(beta);
    const double c = std::cos(beta);

    return {
        beta,
        s, -c,
        c, s
    };
}

inline ThroatLocalVelocity to_throat_local(
    const ThroatVelocityBasis& basis,
    const ImpellerRelativeVelocity& velocity)
{
    if (!std::isfinite(velocity.meridional) ||
        !std::isfinite(velocity.circumferential)) {
        throw std::invalid_argument(
            "Impeller relative velocity must be finite");
    }

    const double Wn =
        basis.normal_meridional * velocity.meridional +
        basis.normal_circumferential *
            velocity.circumferential;

    const double Wt =
        basis.tangent_meridional * velocity.meridional +
        basis.tangent_circumferential *
            velocity.circumferential;

    return {Wn, Wt};
}

inline ImpellerRelativeVelocity from_throat_local(
    const ThroatVelocityBasis& basis,
    const ThroatLocalVelocity& velocity)
{
    if (!std::isfinite(velocity.normal) ||
        !std::isfinite(velocity.tangential)) {
        throw std::invalid_argument(
            "Throat-local velocity must be finite");
    }

    const double Wm =
        basis.normal_meridional * velocity.normal +
        basis.tangent_meridional * velocity.tangential;

    const double Wtheta =
        basis.normal_circumferential * velocity.normal +
        basis.tangent_circumferential *
            velocity.tangential;

    return {Wm, Wtheta};
}

inline ImpellerRelativeVelocity blade_aligned_velocity(
    double blade_angle_deg,
    double relative_speed)
{
    if (!std::isfinite(relative_speed) ||
        relative_speed < 0.0) {
        throw std::invalid_argument(
            "Blade-aligned relative speed must be nonnegative and finite");
    }

    const auto basis =
        make_throat_velocity_basis(blade_angle_deg);

    return from_throat_local(
        basis,
        {relative_speed, 0.0});
}

inline double relative_speed_squared(
    const ImpellerRelativeVelocity& velocity)
{
    return velocity.meridional * velocity.meridional +
           velocity.circumferential *
               velocity.circumferential;
}

inline double relative_speed_squared(
    const ThroatLocalVelocity& velocity)
{
    return velocity.normal * velocity.normal +
           velocity.tangential * velocity.tangential;
}

} // namespace turbocompgpu
