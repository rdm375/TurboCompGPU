#include <turbocompgpu/counting_model.hpp>
#include <turbocompgpu/nonlinear_solver.hpp>

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

struct ProfileCase {
    const char* model_name;
    double entropy_generation;
    double reference_temperature;
    double reference_pressure;
};

constexpr std::array<ProfileCase, 6> cases {{
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

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
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
        300.0,
        500000.0,
        5.0,
        12000.0,
        0.90,
        20.0,
        loss
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

void profile_case(
    const thermomodel::Model& underlying,
    const ProfileCase& test)
{
    CountingModel counted(underlying);

    const auto s = scales();

    const NonlinearUnknowns initial {
        340.0,
        std::log(700000.0 / s.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    counted.reset();

    const auto result =
        solve_nonlinear_compressor(
            counted,
            geometry(),
            operating(test.entropy_generation),
            initial,
            s,
            options);

    const auto counts = counted.counts();

    require(
        result.converged(),
        std::string("Profiling solve failed: ") +
            test.model_name);

    require(
        result.solution.has_value(),
        "Missing profiling solution");

    const auto& solution = *result.solution;

    require(
        std::abs(
            solution.unknowns.outlet_temperature -
            test.reference_temperature) < 1.0e-6,
        "Profiling changed outlet temperature");

    require(
        std::abs(
            solution.outlet_pressure -
            test.reference_pressure) < 1.0e-3,
        "Profiling changed outlet pressure");

    const double calls_per_residual =
        static_cast<double>(counts.total()) /
        result.residual_evaluations;

    std::cout
        << std::left
        << std::setw(15) << test.model_name
        << std::right
        << std::setw(7) << test.entropy_generation
        << std::setw(8) << result.iterations
        << std::setw(11) << result.residual_evaluations
        << std::setw(12) << counts.evaluate_Tp
        << std::setw(12) << counts.evaluate_Trho
        << std::setw(12) << counts.total()
        << std::setw(14) << std::fixed
        << std::setprecision(2) << calls_per_residual
        << '\n';
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

        std::cout
            << std::left
            << std::setw(15) << "Model"
            << std::right
            << std::setw(7) << "Loss"
            << std::setw(8) << "Newton"
            << std::setw(11) << "Residuals"
            << std::setw(12) << "Tp"
            << std::setw(12) << "Trho"
            << std::setw(12) << "Total"
            << std::setw(14) << "Calls/Res"
            << '\n';

        for (const auto& test : cases) {
            const thermomodel::Model& model =
                std::string(test.model_name) == "IdealGas"
                    ? static_cast<const thermomodel::Model&>(
                          ideal_gas)
                    : static_cast<const thermomodel::Model&>(
                          methane);

            profile_case(model, test);
        }

        std::cout
            << "PASS: Thermodynamic evaluation profiling\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
