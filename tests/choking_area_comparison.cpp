#include <turbocompgpu/choking_area_comparison.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double atol,
    double rtol,
    const char* message)
{
    require(
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            atol + rtol * std::max(
                std::abs(actual), std::abs(expected)),
        message);
}

int main()
{
    try {
        thermomodel::IdealGas air;

        const ImpellerSystemGeometry geometry{
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {0.40, 45.0, 0.15}
        };

        const ThroatGeometry throat{
            12, 0.15, 0.0035, 45.0, 0.002, 0.10
        };

        const PassageTurningChokingInputs baseline{
            300.0, 500000.0, 12000.0, 20.0, throat, 0.5
        };

        const auto reference_geometry =
            evaluate_throat_geometry(throat);

        const auto points =
            evaluate_choking_area_comparison(
                air, geometry, baseline,
                {25.0, 45.0, 65.0});

        require(points.size() == 3, "Wrong point count");

        for (const auto& point : points) {
            close(
                point.constant_area,
                reference_geometry.effective_area,
                1e-13, 1e-12,
                "Constant throat area");

            require(
                point.physical_mass_flow > 0.0 &&
                point.constant_area_mass_flow > 0.0,
                "Nonpositive choking capacity");

            require(
                std::abs(point.physical_mass_residual) < 1e-8 &&
                std::abs(point.constant_area_mass_residual) < 1e-8,
                "Outer mass-flow residual");

            require(
                std::abs(point.physical_rothalpy_residual) < 1e-5 &&
                std::abs(point.constant_area_rothalpy_residual) < 1e-5,
                "Rothalpy conservation");
        }

        const auto& reference = points[1];

        close(
            reference.physical_height,
            reference.constant_area_height,
            1e-13, 1e-12,
            "Reference passage height");

        close(
            reference.physical_mass_flow,
            reference.constant_area_mass_flow,
            1e-9, 1e-10,
            "Reference choking mass flow");

        close(
            reference.physical_tangential_velocity,
            reference.constant_area_tangential_velocity,
            1e-9, 1e-10,
            "Reference tangential velocity");

        // The physical-geometry area increases with angle
        // for the selected geometry.
        require(
            points[0].physical_area <
                reference.physical_area &&
            reference.physical_area <
                points[2].physical_area,
            "Unexpected physical-area ordering");

        // At fixed area, the required passage height
        // decreases as blade angle increases.
        require(
            points[0].constant_area_height >
                reference.constant_area_height &&
            reference.constant_area_height >
                points[2].constant_area_height,
            "Unexpected constant-area height ordering");

        std::cout
            << "PASS: choking constant-area comparison\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
