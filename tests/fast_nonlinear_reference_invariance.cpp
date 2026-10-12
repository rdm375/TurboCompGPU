#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/reference_shift_model.hpp>

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

constexpr std::array<double, 3> enthalpy_shifts {{
    -1.0e6, 0.0, 1.0e6
}};

constexpr std::array<double, 3> entropy_shifts {{
    -1.0e4, 0.0, 1.0e4
}};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(
    double actual,
    double expected,
    double relative_tolerance,
    const std::string& message)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <=
            relative_tolerance * scale,
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

struct SolutionSummary {
    double temperature;
    double pressure;
    double euler_work;
    double torque;
    int iterations;
};

SolutionSummary solve_case(
    const thermomodel::Model& model,
    const ReferenceCase& reference)
{
    const auto g = geometry();
    const auto op = operating(reference.entropy_generation);
    const auto s = scales();

    const auto context =
        make_fast_residual_context(model, g, op, s);

    const NonlinearUnknowns initial {
        340.0,
        std::log(700000.0 / s.pressure)
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto result =
        solve_fast_nonlinear_compressor(
            context, initial, options);

    require(
        result.converged(),
        std::string("Fast Newton failed: ") +
            reference.model_name + ": " + result.message);

    require(
        result.solution.has_value(),
        "Converged solver has no solution");

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
        result.final_residual_norm <=
            options.residual_tolerance,
        "Normalized residual tolerance not satisfied");

    // Independent full-system reconstruction.
    const auto validated =
        validate_fast_solution(context, solution);

    const auto& system = validated.system;
    const auto& loss = validated.loss;

    const double mass_scale = 5.0;

    const double torque_scale =
        std::max(1.0, std::abs(solution.shaft_torque));

    constexpr double pi =
        3.141592653589793238462643383279502884;

    const double omega =
        2.0 * pi * op.rotational_speed_rpm / 60.0;

    const double power_scale =
        std::max(
            1.0,
            std::abs(solution.shaft_torque * omega));

    const double whirl_scale =
        std::max(1.0, std::abs(solution.outlet_whirl));

    const double entropy_scale =
        std::max(1.0, std::abs(op.entropy_generation));

    require(
        std::abs(system.conservation.mass) <=
            1.0e-12 * mass_scale,
        "Independent mass conservation failed");

    require(
        std::abs(system.conservation.angular_momentum) <=
            1.0e-12 * torque_scale,
        "Independent angular momentum conservation failed");

    require(
        std::abs(system.conservation.energy) <=
            1.0e-10 * power_scale,
        "Independent energy conservation failed");

    require(
        std::abs(system.slip_residual) <=
            1.0e-12 * whirl_scale,
        "Independent slip closure failed");

    require(
        std::abs(loss.entropy_residual) <=
            1.0e-10 * entropy_scale,
        "Independent entropy closure failed");

    require(
        loss.actual_entropy_increase >= -1.0e-10,
        "Second-law inequality failed");

    return {
        solution.unknowns.outlet_temperature,
        solution.outlet_pressure,
        solution.euler_specific_work,
        solution.shaft_torque,
        result.iterations
    };
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane(caloric);

        int total_solves = 0;
        int changed_iteration_counts = 0;

        double maximum_temperature_difference = 0.0;
        double maximum_pressure_difference = 0.0;
        double maximum_work_difference = 0.0;
        double maximum_torque_difference = 0.0;

        std::cout << std::setprecision(12);

        for (const auto& reference : references) {

            const thermomodel::Model& base =
                std::string(reference.model_name) == "IdealGas"
                    ? static_cast<const thermomodel::Model&>(
                          ideal)
                    : static_cast<const thermomodel::Model&>(
                          methane);

            // Unshifted baseline is computed independently.
            const auto baseline =
                solve_case(base, reference);

            for (const double dh : enthalpy_shifts) {
                for (const double ds : entropy_shifts) {

                    ReferenceShiftModel shifted(base, dh, ds);

                    const auto result =
                        solve_case(shifted, reference);

                    const double dT =
                        std::abs(
                            result.temperature -
                            baseline.temperature);

                    const double dp =
                        std::abs(
                            result.pressure -
                            baseline.pressure);

                    const double dw =
                        std::abs(
                            result.euler_work -
                            baseline.euler_work);

                    const double dtorque =
                        std::abs(
                            result.torque -
                            baseline.torque);

                    maximum_temperature_difference =
                        std::max(
                            maximum_temperature_difference,
                            dT);

                    maximum_pressure_difference =
                        std::max(
                            maximum_pressure_difference,
                            dp);

                    maximum_work_difference =
                        std::max(
                            maximum_work_difference,
                            dw);

                    maximum_torque_difference =
                        std::max(
                            maximum_torque_difference,
                            dtorque);

                    if (result.iterations !=
                        baseline.iterations) {
                        ++changed_iteration_counts;
                    }

                    // These compare physical solutions.
                    // They do not demand bitwise equality.
                    near(
                        result.temperature,
                        baseline.temperature,
                        1.0e-8,
                        "Temperature changed under reference shift");

                    near(
                        result.pressure,
                        baseline.pressure,
                        1.0e-8,
                        "Pressure changed under reference shift");

                    near(
                        result.euler_work,
                        baseline.euler_work,
                        1.0e-8,
                        "Euler work changed under reference shift");

                    near(
                        result.torque,
                        baseline.torque,
                        1.0e-8,
                        "Torque changed under reference shift");

                    std::cout
                        << reference.model_name
                        << " loss=" << reference.entropy_generation
                        << " dh=" << dh
                        << " ds=" << ds
                        << " iterations=" << result.iterations
                        << " dT=" << dT
                        << " dp=" << dp
                        << " dw=" << dw
                        << " dtorque=" << dtorque
                        << '\n';

                    ++total_solves;
                }
            }
        }

        require(
            total_solves == 54,
            "Incorrect reference-shift test count");

        std::cout
            << "\nReference-invariance summary\n"
            << "  shifted solves          = "
            << total_solves << '\n'
            << "  unshifted baselines     = "
            << references.size() << '\n'
            << "  changed iteration count = "
            << changed_iteration_counts << '\n'
            << "  maximum |dT| (K)        = "
            << maximum_temperature_difference << '\n'
            << "  maximum |dp| (Pa)       = "
            << maximum_pressure_difference << '\n'
            << "  maximum |dw| (J/kg)     = "
            << maximum_work_difference << '\n'
            << "  maximum |dtorque| (Nm)  = "
            << maximum_torque_difference << '\n'
            << "PASS: Fast nonlinear reference invariance\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
