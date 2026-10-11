#include <turbocompgpu/nonlinear_solver.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

struct ReferenceCase {
    const char* model_name;
    double entropy_generation;
    double temperature;
    double pressure;
};

struct InitialGuess {
    const char* name;
    double temperature;
    double pressure;
};

constexpr std::array<ReferenceCase, 6> references {{
    {"IdealGas", 0.0,
        329.11606992693578, 691421.5273195625},
    {"IdealGas", 20.0,
        328.99653060808760, 644067.8864092113},
    {"IdealGas", 40.0,
        328.86175551202996, 599859.8380954124},
    {"PengRobinson", 0.0,
        313.41641779670454, 600731.9106110390},
    {"PengRobinson", 20.0,
        313.20996035284452, 576979.5527332224},
    {"PengRobinson", 40.0,
        313.00015991161899, 554118.8983365324}
}};

constexpr std::array<InitialGuess, 3> guesses {{
    {"near", 340.0, 700000.0},
    {"moderate", 300.0, 500000.0},
    {"distant", 380.0, 900000.0}
}};

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void near(
    double actual,
    double expected,
    double tolerance,
    const std::string& message)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <= tolerance * scale,
        message);
}

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

NonlinearOperatingPoint operating(double entropy_generation)
{
    return {
        300.0,
        500000.0,
        5.0,
        12000.0,
        0.90,
        20.0,
        entropy_generation
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

void verify_five_residuals(
    const NonlinearResidualResult& r,
    double entropy_generation)
{
    const auto& system = r.full_system.system;
    const auto& loss = r.full_system.loss;

    const double mass_scale = 5.0;

    const double torque_scale =
        std::max(1.0, std::abs(r.shaft_torque));

    const double omega =
        2.0 *
        3.141592653589793238462643383279502884 *
        12000.0 / 60.0;

    const double power_scale =
        std::max(1.0, std::abs(r.shaft_torque * omega));

    const double whirl_scale =
        std::max(1.0, std::abs(r.outlet_whirl));

    const double entropy_scale =
        std::max(1.0, std::abs(entropy_generation));

    require(
        std::abs(system.conservation.mass) <=
            1.0e-12 * mass_scale,
        "Mass conservation failed");

    require(
        std::abs(system.conservation.angular_momentum) <=
            1.0e-12 * torque_scale,
        "Angular momentum conservation failed");

    require(
        std::abs(system.conservation.energy) <=
            1.0e-10 * power_scale,
        "Energy conservation failed");

    require(
        std::abs(system.slip_residual) <=
            1.0e-12 * whirl_scale,
        "Slip closure failed");

    require(
        std::abs(loss.entropy_residual) <=
            1.0e-10 * entropy_scale,
        "Entropy closure failed");

    // At the reversible limit, allow only the established
    // floating-point admissibility tolerance.
    require(
        loss.actual_entropy_increase >= -1.0e-10,
        "Entropy inequality failed");
}

void test_case(
    const thermomodel::Model& model,
    const ReferenceCase& reference,
    const InitialGuess& guess,
    int& total_iterations,
    int& total_evaluations,
    int& total_rejections)
{
    const auto s = scales();

    const NonlinearUnknowns initial {
        guess.temperature,
        std::log(guess.pressure / s.pressure)
    };

    NonlinearSolverOptions options;

    // Extended conservation validation requires tighter
    // convergence than the production default.
    //
    // With entropy_scale = 100 J/(kg K), this corresponds
    // to a nominal entropy residual of 1e-11 J/(kg K).

    options.residual_tolerance = 1.0e-13;

    const auto result =
        solve_nonlinear_compressor(
            model,
            geometry(),
            operating(reference.entropy_generation),
            initial,
            s,
            options);

    std::cout << std::setprecision(12)
              << reference.model_name
              << " loss=" << reference.entropy_generation
              << " guess=" << guess.name
              << " status="
              << nonlinear_status_name(result.status)
              << " iterations=" << result.iterations
              << " evaluations=" << result.residual_evaluations
              << " rejected=" << result.rejected_steps
              << " norm=" << result.final_residual_norm
              << '\n';

    require(
        result.converged(),
        std::string("Newton failed: ") +
            reference.model_name + "/" + guess.name +
            ": " + result.message);

    require(
        result.solution.has_value(),
        "Converged result has no solution");

    const auto& solution = *result.solution;

    near(
        solution.unknowns.outlet_temperature,
        reference.temperature,
        1.0e-8,
        "Manufactured temperature mismatch");

    near(
        solution.outlet_pressure,
        reference.pressure,
        1.0e-8,
        "Manufactured pressure mismatch");

    require(
        result.final_residual_norm <= 1.0e-10,
        "Normalized convergence criterion failed");

    verify_five_residuals(
        solution,
        reference.entropy_generation);

    total_iterations += result.iterations;
    total_evaluations += result.residual_evaluations;
    total_rejections += result.rejected_steps;
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal_gas;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane(caloric);

        int total_iterations = 0;
        int total_evaluations = 0;
        int total_rejections = 0;
        int total_solves = 0;

        for (const auto& reference : references) {
            const thermomodel::Model& model =
                std::string(reference.model_name) == "IdealGas"
                    ? static_cast<const thermomodel::Model&>(
                          ideal_gas)
                    : static_cast<const thermomodel::Model&>(
                          methane);

            for (const auto& guess : guesses) {
                test_case(
                    model,
                    reference,
                    guess,
                    total_iterations,
                    total_evaluations,
                    total_rejections);

                ++total_solves;
            }
        }

        std::cout
            << "\nExtended Newton validation summary\n"
            << "  solves               = "
            << total_solves << '\n'
            << "  total iterations     = "
            << total_iterations << '\n'
            << "  residual evaluations = "
            << total_evaluations << '\n'
            << "  rejected steps       = "
            << total_rejections << '\n'
            << "  mean iterations      = "
            << static_cast<double>(total_iterations) /
               total_solves << '\n'
            << "  mean evaluations     = "
            << static_cast<double>(total_evaluations) /
               total_solves << '\n'
            << "PASS: Extended nonlinear solver validation\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
