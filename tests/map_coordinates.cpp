#include <turbocompgpu/map_coordinates.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace turbocompgpu;

static void require_close(
    double actual,
    double expected,
    double tolerance = 1.0e-12)
{
    const double scale = std::max(1.0, std::abs(expected));

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        throw std::runtime_error("Numerical comparison failed");
    }
}

int main()
{
    const MapReferenceConditions reference {
        288.15,
        101325.0
    };

    // At reference conditions, corrected and physical
    // coordinates must coincide.
    {
        const auto result = corrected_map_coordinates(
            5.0,
            12000.0,
            {288.15, 101325.0},
            reference
        );

        require_close(result.corrected_mass_flow, 5.0);
        require_close(result.corrected_speed_rpm, 12000.0);
    }

    // Temperature four times reference:
    // corrected flow doubles, corrected speed halves.
    {
        const auto result = corrected_map_coordinates(
            5.0,
            12000.0,
            {4.0 * 288.15, 101325.0},
            reference
        );

        require_close(result.corrected_mass_flow, 10.0);
        require_close(result.corrected_speed_rpm, 6000.0);
    }

    // Pressure twice reference:
    // corrected flow halves; speed is unchanged.
    {
        const auto result = corrected_map_coordinates(
            5.0,
            12000.0,
            {288.15, 2.0 * 101325.0},
            reference
        );

        require_close(result.corrected_mass_flow, 2.5);
        require_close(result.corrected_speed_rpm, 12000.0);
    }

    // Pressure-ratio convention must be preserved.
    {
        const auto result = map_pressure_ratio(
            500000.0,
            750000.0,
            MapPressureConvention::StaticToStatic
        );

        require_close(result.value, 1.5);

        if (result.convention !=
            MapPressureConvention::StaticToStatic) {
            throw std::runtime_error(
                "Pressure convention was not preserved");
        }
    }

    // Invalid inputs must fail explicitly.
    {
        bool rejected = false;

        try {
            (void)corrected_map_coordinates(
                5.0,
                12000.0,
                {-300.0, 500000.0},
                reference
            );
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        if (!rejected) {
            throw std::runtime_error(
                "Invalid inlet temperature was accepted");
        }
    }

    std::cout << "Map-coordinate tests passed\n";
}
