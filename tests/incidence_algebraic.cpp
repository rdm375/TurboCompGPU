#include <turbocompgpu/operating_incidence.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void require_close(double actual, double expected)
{
    const double scale = std::max({
        1.0, std::abs(actual), std::abs(expected)
    });

    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) > 1e-11 * scale) {
        throw std::runtime_error(
            "Angle-based and algebraic incidence losses disagree"
        );
    }
}

} // namespace

int main()
{
    try {
        const ImpellerSystemGeometry geometry {
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {0.40, 45.0, 0.15}
        };

        constexpr std::array<double, 3> densities {
            5.80619, 3.21588, 3.25096
        };

        constexpr std::array<double, 5> speeds {
            9000, 10500, 12000, 13500, 15000
        };

        constexpr std::array<double, 7> flows {
            2, 3, 4, 5, 6, 7, 8
        };

        constexpr std::array<double, 3> coefficients {
            0.0, 0.5, 1.0
        };

        int checked = 0;

        for (double density : densities) {
            const OperatingIncidenceConditions design {
                5.0, 12000.0, density, 20.0
            };

            const double blade_angle =
                calibrate_inlet_blade_angle_deg(
                    geometry, design
                );

            const double blade_rad =
                blade_angle * std::numbers::pi / 180.0;

            for (double speed : speeds) {
                for (double flow : flows) {
                    for (double coefficient : coefficients) {
                        const OperatingIncidenceConditions operating {
                            flow, speed, density, 20.0
                        };

                        const auto result =
                            evaluate_operating_incidence_loss(
                                geometry,
                                operating,
                                {blade_angle, coefficient}
                            );

                        const double transverse =
                            result.inlet_meridional_velocity *
                                std::cos(blade_rad)
                            -
                            (result.inlet_blade_speed -
                             operating.inlet_prewhirl) *
                                std::sin(blade_rad);

                        const double expected =
                            0.5 * coefficient *
                            transverse * transverse;

                        require_close(
                            result.loss.specific_loss,
                            expected
                        );

                        ++checked;
                    }
                }
            }
        }

        if (checked != 315)
            throw std::runtime_error("Unexpected test count");

        std::cout
            << "PASS: 315 algebraic incidence comparisons\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
