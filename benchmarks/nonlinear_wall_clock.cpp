#include <turbocompgpu/fast_nonlinear_solver.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

using Clock = std::chrono::steady_clock;

struct Case {
    const char* model;
    double loss;
    double reference_T;
    double reference_p;
};

constexpr std::array<Case, 6> cases {{
    {"IdealGas",      0.0, 329.11606992693578, 691421.5273195625},
    {"IdealGas",     20.0, 328.99653060808760, 644067.8864092113},
    {"IdealGas",     40.0, 328.86175551202996, 599859.8380954124},
    {"PengRobinson",  0.0, 313.41641779670454, 600731.9106110390},
    {"PengRobinson", 20.0, 313.20996035284452, 576979.5527332224},
    {"PengRobinson", 40.0, 313.00015991161899, 554118.8983365324}
}};

enum class Mode {
    Full,
    Fast,
    FastValidated
};

const char* mode_name(Mode mode)
{
    switch (mode) {
    case Mode::Full:          return "full";
    case Mode::Fast:          return "fast";
    case Mode::FastValidated: return "fast_validated";
    }
    throw std::runtime_error("Unknown benchmark mode");
}

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

NonlinearOperatingPoint operating(double loss)
{
    return {
        300.0, 500000.0, 5.0, 12000.0,
        0.90, 20.0, loss
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

void verify_solution(
    double T,
    double p,
    const Case& test)
{
    require(
        std::isfinite(T) && std::isfinite(p),
        "Non-finite benchmark solution");

    require(
        std::abs(T - test.reference_T) < 1.0e-6,
        "Benchmark temperature differs from reference");

    require(
        std::abs(p - test.reference_p) < 1.0e-3,
        "Benchmark pressure differs from reference");
}

// Volatile accumulation ensures the measured results remain
// observable and cannot simply be discarded by the optimizer.
volatile double checksum = 0.0;

void run_once(
    const thermomodel::Model& model,
    const Case& test,
    Mode mode)
{
    const auto g = geometry();
    const auto op = operating(test.loss);
    const auto s = scales();

    const NonlinearUnknowns initial {
        340.0,
        std::log(700000.0 / s.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    if (mode == Mode::Full) {
        const auto result =
            solve_nonlinear_compressor(
                model, g, op, initial, s, options);

        require(result.converged(), "Full solver failed");
        require(result.solution.has_value(), "Missing full solution");

        const auto& solution = *result.solution;

        verify_solution(
            solution.unknowns.outlet_temperature,
            solution.outlet_pressure,
            test);

        checksum = checksum +
            solution.unknowns.outlet_temperature +
            solution.outlet_pressure * 1.0e-6;

        return;
    }

    // Construct a fresh context on every invocation.
    // Context construction is included in both fast timings.
    const auto context =
        make_fast_residual_context(model, g, op, s);

    const auto result =
        solve_fast_nonlinear_compressor(
            context, initial, options);

    require(result.converged(), "Fast solver failed");
    require(result.solution.has_value(), "Missing fast solution");

    const auto& solution = *result.solution;

    verify_solution(
        solution.unknowns.outlet_temperature,
        solution.outlet_pressure,
        test);

    if (mode == Mode::FastValidated) {
        const auto validated =
            validate_fast_solution(context, solution);

        const auto& system = validated.system;

        require(
            std::abs(system.conservation.mass) < 1.0e-10,
            "Validation: mass conservation");

        require(
            std::abs(system.conservation.angular_momentum) < 1.0e-8,
            "Validation: angular momentum");

        require(
            std::abs(system.conservation.energy) < 1.0e-4,
            "Validation: energy conservation");

        require(
            std::abs(system.slip_residual) < 1.0e-10,
            "Validation: slip closure");

        require(
            std::abs(validated.loss.entropy_residual) < 1.0e-8,
            "Validation: entropy closure");

        require(
            validated.loss.actual_entropy_increase >= -1.0e-10,
            "Validation: entropy inequality");

        checksum = checksum +
            system.conservation.energy;
    }

    checksum = checksum +
        solution.unknowns.outlet_temperature +
        solution.outlet_pressure * 1.0e-6;
}

std::size_t parse_count(const char* value)
{
    const std::string text(value);

    require(!text.empty() && text[0] != '-',
            "Count must be positive");

    std::size_t consumed = 0;
    const auto parsed = std::stoull(text, &consumed);

    require(consumed == text.size() && parsed > 0,
            "Invalid positive count");

    return static_cast<std::size_t>(parsed);
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::string output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-wall-clock/timings.csv";

        const std::size_t repetitions =
            argc > 2 ? parse_count(argv[2]) : 1000;

        const std::size_t warmups =
            argc > 3 ? parse_count(argv[3]) : 50;

        const thermomodel::IdealGas ideal;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        const thermomodel::PengRobinson methane(caloric);

        std::ofstream csv(output);

        require(csv.good(), "Cannot open output CSV");

        csv << "model,loss,mode,repetition,elapsed_ns\n";
        csv << std::setprecision(17);

        constexpr std::array<Mode, 3> modes {{
            Mode::Full,
            Mode::Fast,
            Mode::FastValidated
        }};

        for (const auto& test : cases) {
            const thermomodel::Model& model =
                std::string(test.model) == "IdealGas"
                    ? static_cast<const thermomodel::Model&>(ideal)
                    : static_cast<const thermomodel::Model&>(methane);

            // Warm up each execution path.
            for (std::size_t i = 0; i < warmups; ++i) {
                for (const auto mode : modes)
                    run_once(model, test, mode);
            }

            // Rotate the execution order each repetition.
            // This reduces systematic bias from always
            // timing one implementation first.
            for (std::size_t i = 0; i < repetitions; ++i) {
                for (std::size_t j = 0; j < modes.size(); ++j) {
                    const auto mode = modes[(i + j) % modes.size()];

                    const auto start = Clock::now();

                    run_once(model, test, mode);

                    const auto stop = Clock::now();

                    const auto elapsed =
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(stop - start).count();

                    csv
                        << test.model << ','
                        << test.loss << ','
                        << mode_name(mode) << ','
                        << i << ','
                        << elapsed << '\n';
                }
            }

            std::cout
                << "Completed "
                << test.model
                << " loss=" << test.loss
                << " repetitions=" << repetitions
                << '\n';
        }

        csv.flush();
        require(csv.good(), "CSV write failed");

        std::cout
            << "CSV: " << output << '\n'
            << "Checksum: " << checksum << '\n'
            << "PASS: Wall-clock benchmark completed\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
