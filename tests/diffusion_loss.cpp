#include <turbocompgpu/diffusion_loss.hpp>
#include <turbocompgpu/operating_diffusion.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void check(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(double actual, double expected)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    check(
        std::isfinite(actual) &&
        std::abs(actual - expected) <= 1e-12 * scale,
        "Numerical mismatch"
    );
}

} // namespace

int main()
{
    try {
        constexpr DiffusionLossParameters params {
            1.0, 1.5
        };

        // No diffusion.
        {
            const auto result =
                evaluate_diffusion_loss(100.0, 100.0, params);

            close(result.diffusion_ratio, 1.0);
            close(result.specific_loss, 0.0);
        }

        // Relative-flow acceleration.
        {
            const auto result =
                evaluate_diffusion_loss(100.0, 125.0, params);

            close(result.diffusion_ratio, 0.8);
            close(result.specific_loss, 0.0);
        }

        // Diffusion below the threshold.
        {
            const auto result =
                evaluate_diffusion_loss(100.0, 80.0, params);

            close(result.diffusion_ratio, 1.25);
            close(result.specific_loss, 0.0);
        }

        // Exactly at the threshold.
        {
            const auto result =
                evaluate_diffusion_loss(
                    150.0, 100.0, params
                );

            close(result.diffusion_ratio, 1.5);
            close(result.specific_loss, 0.0);
        }

        // Excess diffusion:
        // W1 = 200, W2 = 100, ratio = 2.
        // Loss = 0.5 * 200^2 * (2 - 1.5)^2
        //      = 5000 J/kg.
        {
            const auto result =
                evaluate_diffusion_loss(
                    200.0, 100.0, params
                );

            close(result.diffusion_ratio, 2.0);
            close(result.excess_diffusion, 0.5);
            close(result.specific_loss, 5000.0);
        }

        // Zero coefficient must recover zero loss.
        {
            const auto result =
                evaluate_diffusion_loss(
                    200.0, 100.0, {0.0, 1.5}
                );

            close(result.specific_loss, 0.0);
        }

        // Quadratic coefficient-independent velocity scaling.
        {
            const auto original =
                evaluate_diffusion_loss(
                    200.0, 100.0, params
                );

            const auto scaled =
                evaluate_diffusion_loss(
                    400.0, 200.0, params
                );

            close(
                scaled.specific_loss,
                4.0 * original.specific_loss
            );
        }

        // Operating-point geometry:
        // Construct velocities so W1 = 200 and W2 = 100.
        {
            const ImpellerSystemGeometry geometry {
                {0.10, 0.20, 0.05},
                {0.40, 0.025, 0.10},
                {0.40, 45.0, 0.15}
            };

            constexpr double rpm = 12000.0;

            const double omega =
                2.0 * std::numbers::pi *
                rpm / 60.0;

            const double U1 =
                0.5 * geometry.impeller.inlet_diameter *
                omega;

            const double U2 =
                0.5 * geometry.impeller.exit_diameter *
                omega;

            const OperatingDiffusionConditions operating {
                200.0,
                100.0,
                U1,
                U2,
                rpm
            };

            const auto result =
                evaluate_operating_diffusion_loss(
                    geometry, operating, params
                );

            close(
                result.loss.inlet_relative_velocity,
                200.0
            );

            close(
                result.loss.outlet_relative_velocity,
                100.0
            );

            close(result.loss.specific_loss, 5000.0);
        }

        // Invalid input rejection.
        {
            bool rejected = false;

            try {
                (void)evaluate_diffusion_loss(
                    100.0, 0.0, params
                );
            }
            catch (const std::invalid_argument&) {
                rejected = true;
            }

            check(
                rejected,
                "Zero outlet relative velocity not rejected"
            );
        }

        std::cout
            << "PASS: diffusion-loss kernel and geometry\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
