#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/nonlinear_solver.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
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

struct Case {
    const char* model;
    double loss;
    double reference_T;
    double reference_p;
};

constexpr std::array<Case, 6> cases {{
    {"IdealGas", 0.0, 329.11606992693578, 691421.5273195625},
    {"IdealGas", 20.0, 328.99653060808760, 644067.8864092113},
    {"IdealGas", 40.0, 328.86175551202996, 599859.8380954124},
    {"PengRobinson", 0.0, 313.41641779670454, 600731.9106110390},
    {"PengRobinson", 20.0, 313.20996035284452, 576979.5527332224},
    {"PengRobinson", 40.0, 313.00015991161899, 554118.8983365324}
}};

constexpr std::array<double, 8> temperatures {{
    260, 280, 300, 320, 340, 360, 380, 400
}};

constexpr std::array<double, 8> pressures_kPa {{
    300, 400, 500, 600, 700, 800, 900, 1000
}};

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

struct Statistics {
    int total = 0;
    int converged = 0;
    int recovered = 0;
    int wrong_branch = 0;
    int failed = 0;
    int iterations = 0;
    int evaluations = 0;
    int rejected = 0;
};

template <typename Result>
void write_result(
    std::ofstream& csv,
    const Case& test,
    const char* solver,
    double initial_T,
    double initial_p,
    const Result& result,
    Statistics& stats)
{
    ++stats.total;

    stats.iterations += result.iterations;
    stats.evaluations += result.residual_evaluations;
    stats.rejected += result.rejected_steps;

    const bool converged = result.converged();

    const double final_T =
        result.solution
            ? result.solution->unknowns.outlet_temperature
            : std::numeric_limits<double>::quiet_NaN();

    const double final_p =
        result.solution
            ? result.solution->outlet_pressure
            : std::numeric_limits<double>::quiet_NaN();

    const double dT = std::abs(final_T - test.reference_T);
    const double dp = std::abs(final_p - test.reference_p);

    // Classification tolerances are not Newton tolerances.
    // They only identify recovery of the manufactured root.
    const bool recovered =
        converged &&
        result.solution.has_value() &&
        dT <= 1.0e-6 &&
        dp <= 1.0e-3;

    if (converged) {
        ++stats.converged;

        if (recovered)
            ++stats.recovered;
        else
            ++stats.wrong_branch;
    }
    else {
        ++stats.failed;
    }

    csv
        << test.model << ','
        << test.loss << ','
        << solver << ','
        << initial_T << ','
        << initial_p << ','
        << nonlinear_status_name(result.status) << ','
        << (converged ? 1 : 0) << ','
        << (recovered ? 1 : 0) << ','
        << result.iterations << ','
        << result.residual_evaluations << ','
        << result.jacobian_evaluations << ','
        << result.rejected_steps << ','
        << result.final_residual_norm << ','
        << final_T << ','
        << final_p << ','
        << dT << ','
        << dp << '\n';
}

void print_summary(
    const char* model,
    const char* solver,
    const Statistics& s)
{
    std::cout
        << std::left << std::setw(16) << model
        << std::setw(8) << solver
        << std::right
        << std::setw(8) << s.total
        << std::setw(10) << s.converged
        << std::setw(10) << s.recovered
        << std::setw(10) << s.wrong_branch
        << std::setw(10) << s.failed
        << std::setw(12) << s.iterations
        << std::setw(12) << s.evaluations
        << std::setw(12) << s.rejected
        << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-convergence-basin/basin.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open output CSV");

        csv << std::setprecision(17);

        csv
            << "model,loss,solver,initial_T,initial_p,"
            << "status,converged,recovered,"
            << "iterations,evaluations,jacobians,rejected,"
            << "residual_norm,final_T,final_p,dT,dp\n";

        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        const thermomodel::PengRobinson methane(caloric);

        NonlinearSolverOptions options;
        options.residual_tolerance = 1.0e-13;

        std::array<std::array<Statistics, 2>, 2> summary {};

        for (const auto& test : cases) {

            const int model_index =
                std::string(test.model) == "IdealGas" ? 0 : 1;

            const thermomodel::Model& model =
                model_index == 0
                    ? static_cast<const thermomodel::Model&>(ideal)
                    : static_cast<const thermomodel::Model&>(methane);

            const auto g = geometry();
            const auto op = operating(test.loss);
            const auto s = scales();

            const auto context =
                make_fast_residual_context(model, g, op, s);

            for (const double initial_T : temperatures) {
                for (const double initial_kPa : pressures_kPa) {

                    const double initial_p = initial_kPa * 1000.0;

                    const NonlinearUnknowns initial {
                        initial_T,
                        std::log(initial_p / s.pressure)
                    };

                    const auto full =
                        solve_nonlinear_compressor(
                            model, g, op, initial, s, options);

                    const auto fast =
                        solve_fast_nonlinear_compressor(
                            context, initial, options);

                    write_result(
                        csv, test, "full",
                        initial_T, initial_p, full,
                        summary[model_index][0]);

                    write_result(
                        csv, test, "fast",
                        initial_T, initial_p, fast,
                        summary[model_index][1]);
                }
            }
        }

        csv.close();

        std::cout
            << std::left
            << std::setw(16) << "Model"
            << std::setw(8) << "Solver"
            << std::right
            << std::setw(8) << "Total"
            << std::setw(10) << "Converged"
            << std::setw(10) << "Recovered"
            << std::setw(10) << "OtherRoot"
            << std::setw(10) << "Failed"
            << std::setw(12) << "Iterations"
            << std::setw(12) << "Evaluations"
            << std::setw(12) << "Rejected"
            << '\n';

        print_summary(
            "IdealGas", "full", summary[0][0]);

        print_summary(
            "IdealGas", "fast", summary[0][1]);

        print_summary(
            "PengRobinson", "full", summary[1][0]);

        print_summary(
            "PengRobinson", "fast", summary[1][1]);

        std::cout
            << "\nCSV: " << output << '\n'
            << "PASS: Convergence-basin experiment completed\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
