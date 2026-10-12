#include <turbocompgpu/choking_continuation.hpp>

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
            300.0, 500000.0, 9000.0, 20.0, throat, 0.5
        };

        const auto points =
            evaluate_choking_continuation(
                air,
                geometry,
                baseline,
                ChokingSweepParameter::turning_fraction,
                {0.0, 0.5, 1.0});

        require(points.size() == 3, "Wrong continuation size");

        for (std::size_t i = 0; i < points.size(); ++i) {
            const auto& point = points[i];

            require(
                point.critical_mass_flow > 0.0,
                "Nonpositive critical mass flow");

            require(
                point.throat_meridional_velocity > 0.0,
                "Reverse sonic throat flow");

            require(
                std::abs(point.mass_flow_residual) < 1e-8,
                "Mass-flow residual");

            require(
                std::abs(point.rothalpy_residual) < 1e-5,
                "Rothalpy residual");

            require(
                std::abs(point.continuity_residual) < 1e-8,
                "Continuity residual");

            require(
                std::abs(point.sonic_residual) < 1e-8,
                "Sonic residual");

            require(
                std::abs(point.implicit_denominator) > 0.9,
                "Unexpected near-singular branch");

            if (i > 0) {
                require(
                    std::abs(
                        point.predicted_mass_flow -
                        points[i-1].critical_mass_flow) < 1e-12,
                    "Predictor not previous converged mass flow");
            }
        }

        require(
            std::abs(points.front().predictor_error) < 1e-12,
            "Initial predictor error must be zero");

        require(
            points.front().critical_mass_flow >
                points.back().critical_mass_flow,
            "Expected test-case choking capacity reduction");

        // Independent direct solve at the middle point.
        const auto direct =
            evaluate_passage_turning_choking(
                air, air, geometry, baseline);

        require(
            std::abs(
                direct.choking.coupled.critical_mass_flow -
                points[1].critical_mass_flow) < 1e-9,
            "Continuation/direct solve disagreement");

        std::cout
            << "PASS: choking continuation diagnostics\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
