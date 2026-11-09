#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace turbocompgpu;

namespace {

struct Point {
    double rpm{};
    double mass_flow{};
    double measured_pr{};
    double T0{};
    double p0{};
    std::string source_run;
    std::string station;
};

std::vector<std::string> split_csv(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}

std::vector<Point> read_points(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open local RADIVER map: " + path.string());
    std::vector<Point> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        const auto f = split_csv(line);
        if (f.size() != 7)
            throw std::runtime_error("Expected 7 columns in local RADIVER CSV: " + line);
        Point point{std::stod(f[0]), std::stod(f[1]), std::stod(f[2]),
                    std::stod(f[3]), std::stod(f[4]), f[5], f[6]};
        if (!(point.rpm > 0.0 && point.mass_flow > 0.0 &&
              point.measured_pr > 1.0 && point.T0 > 0.0 && point.p0 > 0.0))
            throw std::runtime_error("Invalid positive operating data: " + line);
        if (point.station != "impeller_exit_total" &&
            point.station != "vaneless_exit_total")
            throw std::runtime_error(
                "station must be impeller_exit_total or vaneless_exit_total: " + point.station);
        points.push_back(std::move(point));
    }
    if (points.empty()) throw std::runtime_error("Local RADIVER map contains no data rows");
    return points;
}

ImpellerSystemGeometry radiver_geometry()
{
    // Public geometry summary: Jung et al., Table 1. Blade blockage is set to
    // zero because the public summary does not specify blade thickness.
    constexpr double hub_radius = 0.0300;
    constexpr double shroud_radius = 0.0729;
    constexpr double exit_radius = 0.1350;
    constexpr double exit_width = 0.0110;
    constexpr double reported_backsweep_deg = 38.0;
    const double inlet_mean_radius =
        std::sqrt(0.5 * (hub_radius*hub_radius + shroud_radius*shroud_radius));
    return {
        {2.0*hub_radius, 2.0*shroud_radius, 0.0},
        {2.0*exit_radius, exit_width, 0.0},
        // TurboCompGPU's convention is complementary to the reported angle.
        {2.0*exit_radius, 90.0 - reported_backsweep_deg,
         2.0*inlet_mean_radius}
    };
}

struct InletState {
    double T{};
    double p{};
    double velocity{};
};

InletState inlet_from_total(const thermomodel::IdealGas& air,
                            const Point& point,
                            const ImpellerSystemGeometry& geometry)
{
    const double gamma = air.parameters().cp /
        (air.parameters().cp - air.parameters().R);
    const double area = inlet_annulus_area(geometry.inlet);
    auto mass_flow_at_mach = [&](double mach) {
        const double factor = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
        return point.p0 * area / std::sqrt(point.T0) *
            std::sqrt(gamma / air.parameters().R) * mach *
            std::pow(factor, -(gamma + 1.0) / (2.0 * (gamma - 1.0)));
    };
    if (point.mass_flow >= mass_flow_at_mach(1.0))
        throw std::runtime_error("RADIVER inlet flow is at or above the assumed annulus choke limit");
    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (mass_flow_at_mach(mid) < point.mass_flow) lo = mid;
        else hi = mid;
    }
    const double mach = 0.5 * (lo + hi);
    const double factor = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    return {
        point.T0 / factor,
        point.p0 * std::pow(1.0 / factor, gamma / (gamma - 1.0)),
        mach * std::sqrt(gamma * air.parameters().R * point.T0 / factor)
    };
}

std::pair<double, double> stage_total_pressure_ratios(
    const thermomodel::IdealGas& air,
    const ImpellerSystemGeometry& geometry,
    const OperatingPointRequest& request,
    const CompleteStageResult& result)
{
    const auto& stage = *result.stage;
    const auto& station3 = *stage.station3;
    const auto state1 = evaluate(air, request.operating.inlet_temperature,
                                 request.operating.inlet_pressure);
    const double area1 = inlet_annulus_area(geometry.inlet);
    const double vm1 = request.operating.mass_flow / (state1.rho * area1);
    const FlowStation flow1{state1, {vm1, request.operating.inlet_prewhirl}};
    const FlowStation flow2{stage.station2_static,
        {stage.station2_meridional_velocity, stage.station2_tangential_velocity}};
    const FlowStation flow3{station3.state,
        {station3.meridional_velocity, station3.tangential_velocity}};
    const double p01 = evaluate_stagnation_state(air, flow1).stagnation_state.p;
    const double p02 = evaluate_stagnation_state(air, flow2).stagnation_state.p;
    const double p03 = evaluate_stagnation_state(air, flow3).stagnation_state.p;
    return {p02 / p01, p03 / p01};
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc != 3)
            throw std::runtime_error(
                "Usage: turbocompgpu_validate_radiver_local_map LOCAL_MAP.csv OUTPUT.csv");
        const auto points = read_points(argv[1]);
        const auto geometry = radiver_geometry();
        const thermomodel::IdealGas air;
        const NonlinearResidualScales scales{50000.0, 100.0, 300000.0};
        NonlinearSolverOptions solver_options;
        solver_options.residual_tolerance = 1.0e-12;

        const std::filesystem::path output(argv[2]);
        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("Cannot write output CSV: " + output.string());
        csv << std::setprecision(12)
            << "source_run,station,rpm,mass_flow_kg_s,measured_PR,predicted_PR,error_percent,"
               "inlet_static_T_K,inlet_static_p_Pa,slip_factor,valid,message\n";

        int valid = 0;
        double abs_sum = 0.0, max_abs = 0.0;
        for (const auto& point : points) {
            const auto inlet = inlet_from_total(air, point, geometry);
            constexpr double beta_reported_deg = 38.0;
            constexpr double blades = 15.0;
            const double slip = 1.0 -
                std::sqrt(std::cos(beta_reported_deg * std::numbers::pi / 180.0)) /
                std::pow(blades, 0.7);
            constexpr double pressure_scale = 300000.0;
            OperatingPointRequest request{
                {inlet.T, inlet.p, point.mass_flow, point.rpm, slip, 0.0, 0.0},
                {inlet.T * 1.15, std::log(1.5 * inlet.p / pressure_scale)}
            };
            const auto stage = evaluate_complete_stage(
                air, air, geometry,
                {geometry.exit.diameter/2.0,
                 2.0 * 0.345,
                 geometry.exit.width, 0.0},
                request, scales, solver_options);
            double predicted = std::numeric_limits<double>::quiet_NaN();
            double error = std::numeric_limits<double>::quiet_NaN();
            std::string message = stage.message;
            bool point_valid = false;
            if (stage.complete()) {
                const auto [p2_ratio, p3_ratio] =
                    stage_total_pressure_ratios(air, geometry, request, stage);
                predicted = point.station == "impeller_exit_total" ? p2_ratio : p3_ratio;
                error = 100.0 * (predicted / point.measured_pr - 1.0);
                point_valid = std::isfinite(error);
            }
            for (char& c : message)
                if (c == ',' || c == '\n' || c == '\r') c = ';';
            csv << point.source_run << ',' << point.station << ',' << point.rpm << ','
                << point.mass_flow << ',' << point.measured_pr << ',' << predicted << ','
                << error << ',' << inlet.T << ',' << inlet.p << ',' << slip << ','
                << point_valid << ',' << message << '\n';
            if (point_valid) {
                ++valid;
                abs_sum += std::abs(error);
                max_abs = std::max(max_abs, std::abs(error));
            }
            std::cout << point.source_run << " " << point.rpm << " rpm, "
                      << point.mass_flow << " kg/s: "
                      << (point_valid ? "PR error=" + std::to_string(error) + "%" :
                          "no valid stage solution: " + message) << '\n';
        }
        std::cout << "RADIVER local screen: " << valid << '/' << points.size()
                  << " valid; MAPE=" << (valid ? abs_sum / valid : 0.0)
                  << "%; max_abs=" << max_abs << "%\nCSV: " << output << '\n';
        // A screen reports solver validity only; it has no accuracy threshold.
        return valid == static_cast<int>(points.size()) ? 0 : 2;
    }
    catch (const std::exception& e) {
        std::cerr << "RADIVER local screen failed: " << e.what() << '\n';
        return 1;
    }
}
