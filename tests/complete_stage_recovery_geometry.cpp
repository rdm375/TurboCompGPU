#include <turbocompgpu/operating_point.hpp>
#include <turbocompgpu/station.hpp>
#include <turbocompgpu/vaneless_diffuser.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
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
    double maximum_velocity_error = 0.0;
    double maximum_mass_error = 0.0;
    double maximum_energy_error = 0.0;
    double maximum_area_error = 0.0;
    int checked = 0;
};

template<class Model>
void validate_row(
    const Model& model,
    const Row& row,
    Diagnostics& diagnostics)
{
    const auto g = geometry();
    const auto op = request(row.mass, row.rpm);

    const auto impeller = evaluate_operating_point(
        model, g, op, scales(), options());

    if (!(impeller.converged &&
          impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::runtime_error("Impeller validation failed");
    }

    const double area2 = effective_flow_area(
        impeller_exit_area(g.exit), g.exit.blockage);

    const double area3 =
        2.0 * std::numbers::pi *
        row.r3 * row.width3 *
        (1.0 - row.blockage3);

    const double area_error =
        std::abs(area3 - row.area3) / area3;

    const double r2 = g.exit.diameter / 2.0;

    if (std::abs(row.r2 - r2) > 1.0e-12)
        throw std::runtime_error("Inlet radius mismatch");

    if (!(area3 > 0.0 && area3 < area2 && row.r3 > r2))
        throw std::runtime_error("Invalid analytical geometry domain");

    const double alpha = area2 / area3;
    const double lambda = row.r3 / r2;

    const double numerator = alpha * alpha - 1.0;
    const double denominator = 1.0 - 1.0 / (lambda * lambda);

    if (!(numerator > 0.0 && denominator > 0.0))
        throw std::runtime_error("Invalid kinematic square root");

    const double K = std::sqrt(numerator / denominator);

    const double sigma = op.operating.slip_factor;

    const double omega =
        2.0 * std::numbers::pi * row.rpm / 60.0;

    const double beta =
        g.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double vm_analytic =
        sigma * omega * r2 /
        (K + sigma / std::tan(beta));

    const double rho2 = impeller.outlet_density;

    if (!(std::isfinite(rho2) && rho2 > 0.0))
        throw std::runtime_error("Invalid outlet density");

    const double vm_numerical =
        row.mass / (rho2 * area2);

    const double whirl2 =
        sigma * (
            omega * r2 -
            vm_numerical / std::tan(beta)
        );

    const double vm3 =
        vm_numerical * area2 / area3;

    const double whirl3 =
        whirl2 * r2 / row.r3;

    const double energy_difference =
        vm_numerical * vm_numerical +
        whirl2 * whirl2 -
        vm3 * vm3 -
        whirl3 * whirl3;

    const double energy_scale =
        std::max(
            1.0,
            vm_numerical * vm_numerical +
            whirl2 * whirl2
        );

    const double velocity_error =
        std::abs(vm_numerical - vm_analytic) /
        vm_numerical;

    const double mass_analytic =
        rho2 * area2 * vm_analytic;

    const double mass_error =
        std::abs(row.mass - mass_analytic) /
        row.mass;

    const double energy_error =
        std::abs(energy_difference) / energy_scale;

    diagnostics.maximum_velocity_error = std::max(
        diagnostics.maximum_velocity_error,
        velocity_error);

    diagnostics.maximum_mass_error = std::max(
        diagnostics.maximum_mass_error,
        mass_error);

    diagnostics.maximum_energy_error = std::max(
        diagnostics.maximum_energy_error,
        energy_error);

    diagnostics.maximum_area_error = std::max(
        diagnostics.maximum_area_error,
        area_error);

    ++diagnostics.checked;

    constexpr double tolerance = 1.0e-6;

    if (!(velocity_error < tolerance &&
          mass_error < tolerance &&
          energy_error < tolerance &&
          area_error < 1.0e-12)) {

        std::cerr
            << "FAIL model=" << row.model
            << " policy=" << row.policy
            << " radius_ratio=" << row.radius_ratio
            << " rpm=" << row.rpm
            << " velocity_error=" << velocity_error
            << " mass_error=" << mass_error
            << " energy_error=" << energy_error
            << " area_error=" << area_error
            << '\n';

        throw std::runtime_error(
            "Analytical recovery validation failed");
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc != 2)
            throw std::runtime_error(
                "Usage: complete_stage_recovery_geometry <CSV>");

        const auto rows = read_rows(argv[1]);

        if (rows.size() != 120)
            throw std::runtime_error(
                "Expected exactly 120 radius configurations");

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
            << "Analytical recovery validation: PASS\n"
            << "Configurations: " << diagnostics.checked << '\n'
            << "Maximum relative velocity error: "
            << diagnostics.maximum_velocity_error << '\n'
            << "Maximum relative mass error: "
            << diagnostics.maximum_mass_error << '\n'
            << "Maximum normalized energy defect: "
            << diagnostics.maximum_energy_error << '\n'
            << "Maximum relative area error: "
            << diagnostics.maximum_area_error << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
