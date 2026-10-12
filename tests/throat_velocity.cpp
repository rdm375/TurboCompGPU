#include <turbocompgpu/throat_velocity.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

void require(
    bool condition,
    const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double tolerance,
    const std::string& message)
{
    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <=
            tolerance *
            std::max(1.0, std::abs(expected)),
        message);
}

void test_orthonormality()
{
    for (double beta : {15.0, 30.0, 45.0, 60.0, 90.0}) {
        const auto b = make_throat_velocity_basis(beta);

        const double nn =
            b.normal_meridional * b.normal_meridional +
            b.normal_circumferential *
                b.normal_circumferential;

        const double tt =
            b.tangent_meridional * b.tangent_meridional +
            b.tangent_circumferential *
                b.tangent_circumferential;

        const double nt =
            b.normal_meridional * b.tangent_meridional +
            b.normal_circumferential *
                b.tangent_circumferential;

        close(nn, 1.0, 1e-13, "Normal vector norm");
        close(tt, 1.0, 1e-13, "Tangent vector norm");
        close(nt, 0.0, 1e-13, "Basis orthogonality");
    }

    std::cout << "PASS: orthonormal basis\n";
}

void test_round_trip()
{
    for (double beta : {15.0, 30.0, 45.0, 60.0, 90.0}) {
        const auto basis = make_throat_velocity_basis(beta);

        for (const auto original : {
                 ImpellerRelativeVelocity{100.0, -50.0},
                 ImpellerRelativeVelocity{0.0, -200.0},
                 ImpellerRelativeVelocity{175.0, 0.0},
                 ImpellerRelativeVelocity{-20.0, 35.0}
             }) {

            const auto local =
                to_throat_local(basis, original);

            const auto recovered =
                from_throat_local(basis, local);

            close(
                recovered.meridional,
                original.meridional,
                1e-12,
                "Meridional round trip");

            close(
                recovered.circumferential,
                original.circumferential,
                1e-12,
                "Circumferential round trip");
        }
    }

    std::cout << "PASS: forward/inverse transformation\n";
}

void test_energy_invariance()
{
    for (double beta : {15.0, 30.0, 45.0, 60.0, 90.0}) {
        const auto basis = make_throat_velocity_basis(beta);

        const ImpellerRelativeVelocity velocity{
            175.0, -80.0
        };

        const auto local =
            to_throat_local(basis, velocity);

        close(
            relative_speed_squared(local),
            relative_speed_squared(velocity),
            1e-12,
            "Relative kinetic-energy invariance");
    }

    std::cout << "PASS: kinetic-energy invariance\n";
}

void test_blade_aligned()
{
    const double beta = 45.0;
    const double W = 300.0;

    const auto velocity =
        blade_aligned_velocity(beta, W);

    const auto basis =
        make_throat_velocity_basis(beta);

    const auto local =
        to_throat_local(basis, velocity);

    const double expected_component =
        W / std::sqrt(2.0);

    close(
        velocity.meridional,
        expected_component,
        1e-12,
        "Blade-aligned meridional velocity");

    close(
        velocity.circumferential,
        -expected_component,
        1e-12,
        "Blade-aligned circumferential velocity");

    close(
        local.normal,
        W,
        1e-12,
        "Blade-aligned normal velocity");

    close(
        local.tangential,
        0.0,
        1e-12,
        "Blade-aligned tangential velocity");

    std::cout << "PASS: blade-aligned flow\n";
}

void test_normal_sonic_condition()
{
    const double beta = 30.0;
    const double a = 320.0;

    const auto velocity =
        blade_aligned_velocity(beta, a);

    const auto basis =
        make_throat_velocity_basis(beta);

    const auto local =
        to_throat_local(basis, velocity);

    close(
        local.normal,
        a,
        1e-12,
        "Normal sonic velocity");

    close(
        local.tangential,
        0.0,
        1e-12,
        "Sonic blade-aligned tangential velocity");

    close(
        velocity.meridional,
        160.0,
        1e-12,
        "Sonic meridional component");

    close(
        velocity.circumferential,
        -a * std::sqrt(3.0) / 2.0,
        1e-12,
        "Sonic circumferential component");

    std::cout << "PASS: blade-aligned normal sonic condition\n";
}

void test_nonzero_tangential_velocity()
{
    const auto basis =
        make_throat_velocity_basis(45.0);

    const ThroatLocalVelocity original{
        300.0, 40.0
    };

    const auto impeller =
        from_throat_local(basis, original);

    const auto recovered =
        to_throat_local(basis, impeller);

    close(
        recovered.normal,
        300.0,
        1e-12,
        "Normal velocity with incidence");

    close(
        recovered.tangential,
        40.0,
        1e-12,
        "Tangential velocity with incidence");

    close(
        relative_speed_squared(impeller),
        300.0 * 300.0 + 40.0 * 40.0,
        1e-12,
        "Energy with nonzero tangential velocity");

    std::cout
        << "PASS: nonzero throat-tangential velocity\n";
}

void test_invalid_inputs()
{
    bool rejected = false;

    try {
        (void)make_throat_velocity_basis(0.0);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Zero blade angle accepted");

    rejected = false;

    try {
        (void)blade_aligned_velocity(45.0, -1.0);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Negative relative speed accepted");

    std::cout << "PASS: invalid input rejection\n";
}

void print_reference_case()
{
    const double beta = 45.0;
    const double a = 320.0;

    const auto basis =
        make_throat_velocity_basis(beta);

    const auto velocity =
        blade_aligned_velocity(beta, a);

    const auto local =
        to_throat_local(basis, velocity);

    std::cout
        << std::setprecision(12)
        << "Reference blade-aligned sonic case:\n"
        << "  blade angle = " << beta << " deg\n"
        << "  sound speed = " << a << " m/s\n"
        << "  Wm          = "
        << velocity.meridional << " m/s\n"
        << "  Wtheta      = "
        << velocity.circumferential << " m/s\n"
        << "  Wn          = "
        << local.normal << " m/s\n"
        << "  Wtau        = "
        << local.tangential << " m/s\n";
}

} // namespace

int main()
{
    try {
        test_orthonormality();
        test_round_trip();
        test_energy_invariance();
        test_blade_aligned();
        test_normal_sonic_condition();
        test_nonzero_tangential_velocity();
        test_invalid_inputs();
        print_reference_case();

        std::cout
            << "PASS: all throat-velocity tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
