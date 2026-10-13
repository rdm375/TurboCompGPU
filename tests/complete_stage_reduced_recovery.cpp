#include <turbocompgpu/recovery_boundary_solver_optimized.hpp>
#include <turbocompgpu/recovery_boundary_solver.hpp>
#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

namespace {

struct Row {
    std::string model;
    std::string policy;
    double radius_ratio;
    double rpm;
    double mass;
    double r2;
    double r3;
    double width3;
    double blockage3;
    double area3;
};

std::vector<std::string> split_csv(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;

    while (std::getline(stream, field, ','))
        fields.push_back(field);

    return fields;
}

std::vector<Row> read_rows(const std::filesystem::path& path)
{
    std::ifstream input(path);

    if (!input)
        throw std::runtime_error("Cannot open CSV: " + path.string());

    std::string line;
    std::getline(input, line);

    const auto header = split_csv(line);

    auto column = [&](const std::string& name) {
        auto it = std::find(header.begin(), header.end(), name);

        if (it == header.end())
            throw std::runtime_error("Missing CSV column: " + name);

        return static_cast<std::size_t>(it - header.begin());
    };

    const auto c_model = column("model");
    const auto c_policy = column("geometry_policy");
    const auto c_ratio = column("radius_ratio");
    const auto c_rpm = column("rpm");
    const auto c_mass = column("recovery_mass_kg_s");
    const auto c_r2 = column("inlet_radius_m");
    const auto c_r3 = column("outlet_radius_m");
    const auto c_width = column("outlet_width_m");
    const auto c_blockage = column("outlet_blockage");
    const auto c_area = column("outlet_area_m2");
    const auto c_converged = column("recovery_converged");

    std::vector<Row> rows;

    while (std::getline(input, line)) {
        if (line.empty())
            continue;

        const auto fields = split_csv(line);

        if (fields.size() != header.size())
            throw std::runtime_error("Malformed CSV row");

        if (fields[c_converged] != "1")
            throw std::runtime_error("Unconverged recovery root");

        rows.push_back({
            fields[c_model],
            fields[c_policy],
            std::stod(fields[c_ratio]),
            std::stod(fields[c_rpm]),
            std::stod(fields[c_mass]),
            std::stod(fields[c_r2]),
            std::stod(fields[c_r3]),
            std::stod(fields[c_width]),
            std::stod(fields[c_blockage]),
            std::stod(fields[c_area])
        });
    }

    return rows;
}

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

OperatingPointRequest request(double mass, double rpm)
{
    return {
        {300.0, 500000.0, mass, rpm, 0.90, 20.0, 20.0},
        {340.0, std::log(700000.0 / 500000.0)}
    };
}

NonlinearSolverOptions options()
{
    NonlinearSolverOptions result;
    result.residual_tolerance = 1.0e-13;
    return result;
}


struct Diagnostics {
    int checked = 0;
    int impeller_evaluations = 0;
    int impeller_residual_evaluations = 0;
    int optimized_impeller_evaluations = 0;
    int optimized_impeller_residual_evaluations = 0;
    double maximum_optimized_mass_error = 0.0;
    double maximum_optimized_residual = 0.0;
    double maximum_relative_mass_error = 0.0;
    double maximum_relative_scalar_residual = 0.0;
};

template<class Model>
void validate_row(
    const Model& model,
    const Row& row,
    Diagnostics& diagnostics)
{
    const auto g = geometry();

    const VanelessDiffuserGeometry diffuser{
        row.r2,
        row.r3,
        row.width3,
        row.blockage3
    };

    RecoveryBoundaryOptions boundary_options;
    boundary_options.minimum_mass_flow = 0.31;
    boundary_options.maximum_mass_flow = 25.0;
    boundary_options.scan_increment = 0.073;
    boundary_options.relative_mass_tolerance = 1.0e-12;
    boundary_options.relative_residual_tolerance = 1.0e-11;

    const auto reduced = solve_recovery_boundary(
        model,
        g,
        diffuser,
        request(5.0, row.rpm),
        scales(),
        options(),
        boundary_options);

    if (!reduced.converged()) {
        std::cerr
            << "Reduced solver failed: "
            << row.model << " "
            << row.policy << " "
            << row.radius_ratio << " "
            << row.rpm << " "
            << reduced.message << '\n';

        throw std::runtime_error(
            "Reduced recovery solver did not converge");
    }

    OptimizedRecoveryBoundaryOptions optimized_options;
    optimized_options.relative_mass_tolerance = 1.0e-12;
    optimized_options.relative_residual_tolerance = 1.0e-11;

    const auto optimized_result = solve_recovery_boundary_optimized(
        model,
        g,
        diffuser,
        request(5.0, row.rpm),
        scales(),
        options(),
        optimized_options);

    if (!optimized_result.converged()) {
        throw std::runtime_error(
            "Optimized recovery failed: " +
            optimized_result.message);
    }

    const double optimized_mass_error =
        std::abs(optimized_result.mass_flow - row.mass) / row.mass;

    const double optimized_residual =
        std::abs(optimized_result.residual) /
        std::max(1.0, optimized_result.mass_flow);

    diagnostics.maximum_optimized_mass_error = std::max(
        diagnostics.maximum_optimized_mass_error,
        optimized_mass_error);

    diagnostics.maximum_optimized_residual = std::max(
        diagnostics.maximum_optimized_residual,
        optimized_residual);

    diagnostics.optimized_impeller_evaluations +=
        optimized_result.impeller_evaluations;

    diagnostics.optimized_impeller_residual_evaluations +=
        optimized_result.impeller_residual_evaluations;

    if (!(optimized_mass_error < 1.0e-8 &&
          optimized_residual < 1.0e-9)) {
        throw std::runtime_error(
            "Optimized recovery disagrees with reference");
    }

    const double mass_error =
        std::abs(reduced.mass_flow - row.mass) / row.mass;

    if (diagnostics.checked < 10) {
        std::cout
            << std::setprecision(17)
            << "model=" << row.model
            << " policy=" << row.policy
            << " radius=" << row.radius_ratio
            << " rpm=" << row.rpm
            << " reference=" << row.mass
            << " reduced=" << reduced.mass_flow
            << " difference=" << reduced.mass_flow - row.mass
            << " relative_error=" << mass_error
            << '\n';
    }

    const double residual_error =
        std::abs(reduced.residual) /
        std::max(1.0, reduced.mass_flow);

    diagnostics.maximum_relative_mass_error = std::max(
        diagnostics.maximum_relative_mass_error,
        mass_error);

    diagnostics.maximum_relative_scalar_residual = std::max(
        diagnostics.maximum_relative_scalar_residual,
        residual_error);

    diagnostics.impeller_evaluations +=
        reduced.impeller_evaluations;

    diagnostics.impeller_residual_evaluations +=
        reduced.impeller_residual_evaluations;

    ++diagnostics.checked;

    constexpr double mass_tolerance = 1.0e-6;
    constexpr double residual_tolerance = 1.0e-9;

    if (!(mass_error < mass_tolerance &&
          residual_error < residual_tolerance)) {
        std::cerr
            << "FAIL model=" << row.model
            << " policy=" << row.policy
            << " ratio=" << row.radius_ratio
            << " rpm=" << row.rpm
            << " reference=" << row.mass
            << " reduced=" << reduced.mass_flow
            << " mass_error=" << mass_error
            << " residual_error=" << residual_error
            << '\n';

        throw std::runtime_error(
            "Reduced recovery disagrees with reference");
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc != 2)
            throw std::runtime_error(
                "Usage: reduced_recovery_reference <CSV>");

        const auto rows = read_rows(argv[1]);

        if (rows.size() != 120)
            throw std::runtime_error(
                "Expected 120 reference configurations");

        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        Diagnostics diagnostics;

        for (const auto& row : rows) {
            if (row.model == "IdealGas") {
                validate_row(ideal, row, diagnostics);
            } else if (row.model == "PengRobinson") {
                validate_row(pr, row, diagnostics);
            } else {
                throw std::runtime_error(
                    "Unknown thermodynamic model");
            }
        }

        std::cout
            << "Reduced recovery reference validation: PASS\\n"
            << "Configurations: " << diagnostics.checked << '\n'
            << "Maximum relative mass error: "
            << diagnostics.maximum_relative_mass_error << '\n'
            << "Maximum relative scalar residual: "
            << diagnostics.maximum_relative_scalar_residual << '\n'
            << "Total impeller evaluations: "
            << diagnostics.impeller_evaluations << '\n'
            << "Total impeller residual evaluations: "
            << diagnostics.impeller_residual_evaluations << '\n';

        std::cout
            << "\nOptimized maximum relative mass error: "
            << diagnostics.maximum_optimized_mass_error
            << "\nOptimized maximum relative scalar residual: "
            << diagnostics.maximum_optimized_residual
            << "\nOptimized impeller evaluations: "
            << diagnostics.optimized_impeller_evaluations
            << "\nOptimized impeller residual evaluations: "
            << diagnostics.optimized_impeller_residual_evaluations
            << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
