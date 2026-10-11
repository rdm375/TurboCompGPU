#include <turbocompgpu/counting_model.hpp>
#include <turbocompgpu/fast_nonlinear_solver.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

struct Case {
    const char* name;
    double loss;
    double temperature;
    double pressure;
};

constexpr std::array<Case, 6> cases {{
    {"IdealGas", 0.0, 329.11606992693578, 691421.5273195625},
    {"IdealGas", 20.0, 328.99653060808760, 644067.8864092113},
    {"IdealGas", 40.0, 328.86175551202996, 599859.8380954124},
    {"PengRobinson", 0.0, 313.41641779670454, 600731.91061103903},
    {"PengRobinson", 20.0, 313.20996035284452, 576979.55273322237},
    {"PengRobinson", 40.0, 313.00015991161899, 554118.89833653241}
}};

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

void run_case(
    const thermomodel::Model& model,
    const Case& test)
{
    const auto s = scales();
    const auto g = geometry();
    const auto op = operating(test.loss);

    const NonlinearUnknowns initial {
        340.0,
        std::log(700000.0 / s.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    CountingModel full_model(model);

    const auto full =
        solve_nonlinear_compressor(
            full_model, g, op, initial, s, options);

    const auto full_counts = full_model.counts();

    require(full.converged(), "Full solver failed");
    require(full.solution.has_value(), "Missing full solution");

    CountingModel fast_model(model);

    const auto context =
        make_fast_residual_context(
            fast_model, g, op, s);

    const auto fast =
        solve_fast_nonlinear_compressor(
            context, initial, options);

    const auto fast_counts = fast_model.counts();

    require(fast.converged(), "Fast solver failed");
    require(fast.solution.has_value(), "Missing fast solution");

    const auto& a = *full.solution;
    const auto& b = *fast.solution;

    require(
        std::abs(b.unknowns.outlet_temperature -
                 test.temperature) < 1.0e-6,
        "Fast temperature differs from manufactured reference");

    require(
        std::abs(b.outlet_pressure -
                 test.pressure) < 1.0e-3,
        "Fast pressure differs from manufactured reference");

    require(
        std::abs(b.unknowns.outlet_temperature -
                 a.unknowns.outlet_temperature) < 1.0e-6,
        "Fast and full temperatures disagree");

    require(
        std::abs(b.outlet_pressure -
                 a.outlet_pressure) < 1.0e-3,
        "Fast and full pressures disagree");

    require(
        std::abs(b.energy_residual) < 1.0e-7,
        "Fast energy residual too large");

    require(
        std::abs(b.entropy_residual) < 1.0e-10,
        "Fast entropy residual too large");

    // Independent final five-residual validation.
    //
    // These calls are deliberately excluded from the
    // Newton-loop evaluation count.

    const auto validated =
        validate_fast_solution(context, b);

    const auto& system = validated.system;

    require(
        std::abs(system.conservation.mass) < 1.0e-10,
        "Mass conservation failed");

    require(
        std::abs(system.conservation.angular_momentum) < 1.0e-8,
        "Angular momentum conservation failed");

    require(
        std::abs(system.conservation.energy) < 1.0e-4,
        "Energy conservation failed");

    require(
        std::abs(system.slip_residual) < 1.0e-10,
        "Slip closure failed");

    require(
        std::abs(validated.loss.entropy_residual) < 1.0e-8,
        "Entropy closure failed");

    require(
        validated.loss.actual_entropy_increase >= -1.0e-10,
        "Entropy inequality failed");

    // The fast context requires one inlet evaluation.
    // Each Newton residual requires one outlet evaluation.
    //
    // All six baseline cases use valid finite-difference
    // and line-search trial states.

    require(
        fast_counts.total() ==
            static_cast<std::size_t>(
                1 + fast.residual_evaluations),
        "Unexpected fast thermodynamic evaluation count");

    const double reduction =
        static_cast<double>(full_counts.total()) /
        static_cast<double>(fast_counts.total());

    std::cout
        << std::left << std::setw(15) << test.name
        << std::right << std::setw(7) << test.loss
        << std::setw(10) << full.iterations
        << std::setw(10) << fast.iterations
        << std::setw(12) << full_counts.total()
        << std::setw(12) << fast_counts.total()
        << std::setw(12) << std::fixed
        << std::setprecision(2) << reduction
        << '\n';
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        const thermomodel::PengRobinson methane(caloric);

        std::cout
            << std::left << std::setw(15) << "Model"
            << std::right << std::setw(7) << "Loss"
            << std::setw(10) << "FullIter"
            << std::setw(10) << "FastIter"
            << std::setw(12) << "FullCalls"
            << std::setw(12) << "FastCalls"
            << std::setw(12) << "Reduction"
            << '\n';

        for (const auto& test : cases) {
            const thermomodel::Model& model =
                std::string(test.name) == "IdealGas"
                    ? static_cast<const thermomodel::Model&>(ideal)
                    : static_cast<const thermomodel::Model&>(methane);

            run_case(model, test);
        }

        std::cout
            << "PASS: Fast nonlinear solver validation\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
