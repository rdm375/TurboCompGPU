#include <turbocompgpu/throat_geometry.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

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

using turbocompgpu::ThroatGeometry;
using turbocompgpu::evaluate_throat_geometry;

void test_zero_thickness_limit()
{
    const ThroatGeometry geometry{
        12, 0.15, 0.025, 45.0, 0.0, 0.0
    };

    const auto result =
        evaluate_throat_geometry(geometry);

    const double expected =
        2.0 * std::acos(-1.0) *
        geometry.radius *
        geometry.passage_height *
        std::sin(std::acos(-1.0) / 4.0);

    close(
        result.geometric_area,
        expected,
        1e-12,
        "Zero-thickness annular projection");

    close(
        result.effective_area,
        expected,
        1e-12,
        "Zero-thickness effective area");

    std::cout
        << "PASS: zero-thickness geometric limit\n";
}

void test_blade_count_invariance()
{
    const ThroatGeometry g1{
        8, 0.15, 0.025, 45.0, 0.0, 0.0
    };

    const ThroatGeometry g2{
        16, 0.15, 0.025, 45.0, 0.0, 0.0
    };

    const auto a = evaluate_throat_geometry(g1);
    const auto b = evaluate_throat_geometry(g2);

    close(
        a.geometric_area,
        b.geometric_area,
        1e-12,
        "Zero-thickness blade-count invariance");

    std::cout
        << "PASS: zero-thickness blade-count invariance\n";
}

void test_thickness_scaling()
{
    const ThroatGeometry geometry{
        12, 0.15, 0.025, 45.0, 0.002, 0.0
    };

    const auto result =
        evaluate_throat_geometry(geometry);

    const double expected =
        2.0 * std::acos(-1.0) *
        geometry.radius *
        geometry.passage_height *
        std::sin(std::acos(-1.0) / 4.0)
        -
        geometry.blade_count *
        geometry.passage_height *
        geometry.blade_thickness;

    close(
        result.geometric_area,
        expected,
        1e-12,
        "Blade-thickness area scaling");

    std::cout
        << "PASS: blade-thickness scaling\n";
}

void test_blockage_scaling()
{
    const ThroatGeometry geometry{
        12, 0.15, 0.025, 45.0, 0.002, 0.15
    };

    const auto result =
        evaluate_throat_geometry(geometry);

    close(
        result.effective_area,
        0.85 * result.geometric_area,
        1e-12,
        "Blockage scaling");

    std::cout
        << "PASS: blockage scaling\n";
}

void test_radius_scaling()
{
    const ThroatGeometry g1{
        12, 0.15, 0.025, 45.0, 0.0, 0.0
    };

    const ThroatGeometry g2{
        12, 0.30, 0.025, 45.0, 0.0, 0.0
    };

    const auto a = evaluate_throat_geometry(g1);
    const auto b = evaluate_throat_geometry(g2);

    close(
        b.geometric_area,
        2.0 * a.geometric_area,
        1e-12,
        "Radius scaling");

    std::cout
        << "PASS: radius scaling\n";
}

void test_invalid_geometry()
{
    auto rejected = [](ThroatGeometry geometry) {
        try {
            (void)evaluate_throat_geometry(geometry);
        }
        catch (const std::invalid_argument&) {
            return true;
        }
        return false;
    };

    const ThroatGeometry base{
        12, 0.15, 0.025, 45.0, 0.002, 0.0
    };

    auto g = base;
    g.blade_count = 1;
    require(rejected(g), "Invalid blade count accepted");

    g = base;
    g.radius = 0.0;
    require(rejected(g), "Zero radius accepted");

    g = base;
    g.passage_height = 0.0;
    require(rejected(g), "Zero passage height accepted");

    g = base;
    g.blade_angle_deg = 0.0;
    require(rejected(g), "Zero blade angle accepted");

    g = base;
    g.blade_thickness = 1.0;
    require(rejected(g), "Closed throat accepted");

    g = base;
    g.blockage = 1.0;
    require(rejected(g), "Complete blockage accepted");

    std::cout
        << "PASS: invalid geometry rejection\n";
}

void test_reference_geometry()
{
    const ThroatGeometry geometry{
        12, 0.15, 0.025, 45.0, 0.002, 0.10
    };

    const auto result =
        evaluate_throat_geometry(geometry);

    std::cout
        << std::setprecision(12)
        << "Reference throat geometry:\n"
        << "  blade count        = "
        << geometry.blade_count << '\n'
        << "  radius             = "
        << geometry.radius << " m\n"
        << "  passage height     = "
        << geometry.passage_height << " m\n"
        << "  blade angle        = "
        << geometry.blade_angle_deg << " deg\n"
        << "  blade thickness    = "
        << geometry.blade_thickness << " m\n"
        << "  pitch              = "
        << result.circumferential_pitch << " m\n"
        << "  projected pitch    = "
        << result.projected_pitch << " m\n"
        << "  open width         = "
        << result.open_width << " m\n"
        << "  geometric area     = "
        << result.geometric_area << " m^2\n"
        << "  effective area     = "
        << result.effective_area << " m^2\n"
        << "  open area fraction = "
        << result.open_area_fraction << '\n';

    require(
        result.effective_area > 0.0,
        "Nonpositive reference throat area");
}

} // namespace

int main()
{
    try {
        test_zero_thickness_limit();
        test_blade_count_invariance();
        test_thickness_scaling();
        test_blockage_scaling();
        test_radius_scaling();
        test_invalid_geometry();
        test_reference_geometry();

        std::cout
            << "PASS: all throat-geometry tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
