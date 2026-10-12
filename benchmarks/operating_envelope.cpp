#include <turbocompgpu/operating_point.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

void write_row(
    std::ofstream& csv,
    const std::string& model,
    const std::string& method,
    const OperatingPointResult& r)
{
    csv << model << ','
        << method << ','
        << r.operating.mass_flow << ','
        << r.operating.rotational_speed_rpm << ','
        << nonlinear_status_name(r.status) << ','
        << r.converged << ','
        << r.independently_validated << ','
        << r.outlet_temperature << ','
        << r.outlet_pressure << ','
        << r.pressure_ratio << ','
        << r.specific_work << ','
        << r.shaft_power << ','
        << r.torque << ','
        << r.entropy_generation << ','
        << r.residual_norm << ','
        << r.validation_residual_norm << ','
        << r.iterations << ','
        << r.evaluations << ','
        << r.jacobians << ','
        << r.rejected_steps << ','
        << r.basic_physics_valid << ','
        << r.inlet_density << ','
        << r.outlet_density << ','
        << (r.operating_limits ==
            OperatingLimitStatus::NotEvaluated
                ? "NotEvaluated" : "Evaluated")
        << '\n';
}

void sweep(
    const thermomodel::Model& model,
    const std::string& model_name,
    std::ofstream& csv)
{
    const auto g = geometry();
    const auto s = scales();

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    constexpr std::array<double, 7> flows {
        2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
    };

    constexpr std::array<double, 5> speeds {
        9000.0, 10500.0, 12000.0, 13500.0, 15000.0
    };

    int converged = 0;
    int validated = 0;
    int total = 0;
    int continuation_better = 0;

    for (double speed : speeds) {
        NonlinearUnknowns previous {
            340.0, std::log(700000.0 / s.pressure)
        };

        bool previous_valid = false;

        for (double flow : flows) {
            NonlinearOperatingPoint operating {
                300.0, 500000.0, flow, speed,
                0.90, 20.0, 20.0
            };

            const NonlinearUnknowns baseline {
                340.0, std::log(700000.0 / s.pressure)
            };

            const auto independent = evaluate_operating_point(
                model, g, {operating, baseline}, s, options);

            const auto continued = evaluate_operating_point(
                model, g,
                {operating, previous_valid ? previous : baseline},
                s, options);

            write_row(csv, model_name, "independent", independent);
            write_row(csv, model_name, "continuation", continued);

            ++total;

            if (independent.converged && continued.converged)
                ++converged;

            if (independent.independently_validated &&
                continued.independently_validated)
                ++validated;

            if (continued.converged &&
                continued.iterations < independent.iterations)
                ++continuation_better;

            if (continued.converged &&
                continued.independently_validated) {
                previous = {
                    continued.outlet_temperature,
                    std::log(continued.outlet_pressure / s.pressure)
                };
                previous_valid = true;
            } else {
                // Do not propagate failed or unvalidated states.
                previous_valid = false;
            }
        }
    }

    std::cout << model_name
              << ": paired=" << total
              << " both_converged=" << converged
              << " both_validated=" << validated
              << " continuation_fewer_iterations="
              << continuation_better << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-operating-envelope/envelope.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open CSV output");

        csv << std::setprecision(17);

        csv << "model,method,mass_flow,rotational_speed_rpm,"
               "status,converged,validated,outlet_temperature,"
               "outlet_pressure,pressure_ratio,specific_work,"
               "shaft_power,torque,entropy_generation,"
               "residual_norm,validation_residual_norm,"
               "iterations,evaluations,jacobians,rejected_steps,"
               "basic_physics_valid,inlet_density,outlet_density,"
               "operating_limits\n";

        const thermomodel::IdealGas ideal;
        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);
        const thermomodel::PengRobinson pr(caloric);

        sweep(ideal, "IdealGas", csv);
        sweep(pr, "PengRobinson", csv);

        std::cout << "CSV: " << output << '\n';
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
