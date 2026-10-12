#include <turbocompgpu/incidence_loss.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace turbocompgpu;

void check_close(double actual, double expected)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (std::abs(actual - expected) > 1e-12 * scale)
        throw std::runtime_error("Incidence-loss mismatch");
}

int main()
{
    try {
        const IncidenceLossParameters params {
            45.0, 1.0
        };

        // Zero incidence must produce zero loss.
        const auto design = evaluate_incidence_loss(
            45.0, 100.0, params
        );

        check_close(design.specific_loss, 0.0);

        // Positive and negative incidence should
        // produce equal losses in this symmetric model.
        const auto positive = evaluate_incidence_loss(
            55.0, 100.0, params
        );

        const auto negative = evaluate_incidence_loss(
            35.0, 100.0, params
        );

        check_close(
            positive.specific_loss,
            negative.specific_loss
        );

        // Doubling velocity quadruples the loss.
        const auto faster = evaluate_incidence_loss(
            55.0, 200.0, params
        );

        check_close(
            faster.specific_loss,
            4.0 * positive.specific_loss
        );

        // Zero coefficient produces zero loss.
        const auto zero = evaluate_incidence_loss(
            55.0, 100.0, {45.0, 0.0}
        );

        check_close(zero.specific_loss, 0.0);

        // Analytical 90-degree case.
        const auto maximum = evaluate_incidence_loss(
            90.0, 100.0, {0.0, 1.0}
        );

        check_close(maximum.specific_loss, 5000.0);

        // Invalid coefficient must be rejected.
        bool rejected = false;

        try {
            (void)evaluate_incidence_loss(
                55.0, 100.0, {45.0, -1.0}
            );
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        if (!rejected)
            throw std::runtime_error(
                "Negative coefficient accepted"
            );

        std::cout
            << "PASS: incidence-loss kernel\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
