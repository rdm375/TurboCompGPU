#include <turbocompgpu/nonlinear_solver.hpp>
#include <turbocompgpu/reference_shift_model.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void near(
    double actual,
    double expected,
    double relative_tolerance,
    const char* message)
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

NonlinearOperatingPoint operating()
{
    return {
        300.0,
        500000.0,
        5.0,
        12000.0,
        0.90,
        20.0,
        20.0
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

NonlinearUnknowns initial_guess()
{
    return {
        340.0,
        std::log(700000.0 / scales().pressure)
    };
}

NonlinearSolverResult run(
    const thermomodel::Model& model,
    const NonlinearUnknowns& guess = initial_guess(),
    const NonlinearSolverOptions& options = {})
{
    return solve_nonlinear_compressor(
        model,
        geometry(),
        operating(),
        guess,
        scales(),
        options);
}

void print_result(
    const char* name,
    const NonlinearSolverResult& r)
{
    std::cout << std::setprecision(17)
              << name << '\n'
              << "  status                = "
              << nonlinear_status_name(r.status) << '\n'
              << "  iterations            = "
              << r.iterations << '\n'
              << "  residual evaluations  = "
              << r.residual_evaluations << '\n'
              << "  Jacobian evaluations  = "
              << r.jacobian_evaluations << '\n'
              << "  rejected steps        = "
              << r.rejected_steps << '\n'
              << "  final residual norm   = "
              << r.final_residual_norm << '\n';

    if (r.solution) {
        std::cout
            << "  T2                    = "
            << r.solution->unknowns.outlet_temperature
            << " K\n"
            << "  p2                    = "
            << r.solution->outlet_pressure
            << " Pa\n"
            << "  energy residual       = "
            << r.solution->energy_residual
            << " J/kg\n"
            << "  entropy residual      = "
            << r.solution->entropy_residual
            << " J/(kg K)\n";
    }
}

void test_manufactured(
    const thermomodel::Model& model,
    double expected_temperature,
    double expected_pressure,
    const char* name)
{
    const auto r = run(model);
    print_result(name, r);

    require(r.converged(), "Newton solver did not converge");
    require(r.solution.has_value(), "Missing converged solution");

    near(
        r.solution->unknowns.outlet_temperature,
        expected_temperature,
        1.0e-8,
        "Incorrect outlet temperature");

    near(
        r.solution->outlet_pressure,
        expected_pressure,
        1.0e-8,
        "Incorrect outlet pressure");

    require(
        r.final_residual_norm <= 1.0e-10,
        "Residual tolerance not satisfied");

    const auto& system = r.solution->full_system.system;

    near(
        system.conservation.mass,
        0.0,
        1.0e-10,
        "Mass conservation failed");

    near(
        system.conservation.angular_momentum,
        0.0,
        1.0e-10,
        "Angular momentum conservation failed");

    near(
        system.slip_residual,
        0.0,
        1.0e-10,
        "Slip closure failed");

    require(
        std::abs(system.conservation.energy) <
            1.0e-4,
        "Energy conservation failed");

    require(
        std::abs(r.solution->full_system.loss.entropy_residual) <
            1.0e-8,
        "Entropy closure failed");

    std::cout << "PASS: " << name << '\n';
}

void test_reference_invariance(
    const thermomodel::Model& model,
    const char* name)
{
    const ReferenceShiftModel shifted(
        model, 1.0e6, 1.0e4);

    const auto a = run(model);
    const auto b = run(shifted);

    require(
        a.converged() && b.converged(),
        "Reference-shifted Newton solve did not converge");

    near(
        b.solution->unknowns.outlet_temperature,
        a.solution->unknowns.outlet_temperature,
        1.0e-8,
        "Temperature depends on reference convention");

    near(
        b.solution->outlet_pressure,
        a.solution->outlet_pressure,
        1.0e-8,
        "Pressure depends on reference convention");

    near(
        b.solution->outlet_whirl,
        a.solution->outlet_whirl,
        1.0e-8,
        "Whirl depends on reference convention");

    near(
        b.solution->shaft_torque,
        a.solution->shaft_torque,
        1.0e-8,
        "Torque depends on reference convention");

    std::cout
        << "PASS: Solver reference invariance: "
        << name << '\n';
}

void test_invalid_initial_guess()
{
    const thermomodel::IdealGas model;

    const auto r = run(
        model,
        {-10.0, 0.0});

    require(
        r.status == NonlinearSolveStatus::InvalidInput,
        "Invalid initial guess was not rejected");

    require(
        !r.solution.has_value(),
        "Invalid guess unexpectedly produced a solution");

    std::cout << "PASS: Invalid initial guess\n";
}

void test_iteration_limit()
{
    const thermomodel::IdealGas model;

    NonlinearSolverOptions options;
    options.maximum_iterations = 0;

    const auto r = run(model, initial_guess(), options);

    require(
        r.status ==
            NonlinearSolveStatus::MaximumIterations,
        "Iteration limit was not reported");

    require(
        r.solution.has_value(),
        "Iteration-limited solve lost its last state");

    require(
        r.iterations == 0,
        "Iteration count incorrect");

    std::cout << "PASS: Iteration limit diagnostics\n";
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

        test_manufactured(
            ideal_gas,
            328.9965306080876,
            644067.88640921132,
            "Ideal-gas Newton solve");

        test_manufactured(
            methane,
            313.20996035284452,
            576979.55273322237,
            "Peng-Robinson Newton solve");

        test_reference_invariance(
            ideal_gas, "Ideal gas");

        test_reference_invariance(
            methane, "Peng-Robinson methane");

        test_invalid_initial_guess();
        test_iteration_limit();

        std::cout
            << "PASS: All nonlinear solver tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
