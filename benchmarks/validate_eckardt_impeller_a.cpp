#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

namespace {

struct MapPoint {
    double rpm;
    double mass_flow;
    double measured_total_pressure_ratio;
};

std::vector<MapPoint> read_map(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open map data: " + path.string());

    std::vector<MapPoint> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        MapPoint point{};
        char comma1 = '\0';
        char comma2 = '\0';
        if (!(row >> point.rpm >> comma1 >> point.mass_flow >> comma2
                  >> point.measured_total_pressure_ratio) ||
            comma1 != ',' || comma2 != ',') {
            throw std::runtime_error("Invalid map row: " + line);
        }
        points.push_back(point);
    }
    if (points.empty())
        throw std::runtime_error("Map data contains no points");
    return points;
}

ImpellerSystemGeometry geometry(bool rotor_b = false)
{
    if (rotor_b) {
        // Eckardt B geometry from Meroni et al., Table 5. The inlet
        // representative diameter is the RMS diameter of the hub/shroud
        // annulus; the model uses a circular-equivalent inlet area.
        constexpr double exit_diameter_m = 0.400;
        constexpr double exit_width_m = 0.026;
        constexpr double blade_count = 20.0;
        constexpr double exit_blade_thickness_m = 0.00108;
        constexpr double exit_blockage = blade_count * exit_blade_thickness_m /
            (2.0 * std::numbers::pi * (exit_diameter_m / 2.0));
        constexpr double inlet_hub_radius_m = 0.0959;
        constexpr double inlet_shroud_radius_m = 0.140;
        return {
            {2.0 * inlet_hub_radius_m, 2.0 * inlet_shroud_radius_m, 0.0},
            {exit_diameter_m, exit_width_m, exit_blockage},
            {exit_diameter_m, 50.0,
             2.0 * std::sqrt(0.5 * (inlet_hub_radius_m * inlet_hub_radius_m +
                                      inlet_shroud_radius_m * inlet_shroud_radius_m))}
        };
    }
    // Eckardt A dimensions transcribed from the public MATLAB model by
    // Guido Fratte. The inlet radii and inlet diameter below use a circular
    // area-equivalent approximation for the annulus; they are not a full
    // reconstruction of the original inlet geometry.
    constexpr double exit_diameter_m = 0.400;
    constexpr double exit_width_m = 0.026;
    constexpr double blade_count = 20.0;
    constexpr double exit_blade_thickness_m = 0.00108;
    constexpr double exit_blockage =
        blade_count * exit_blade_thickness_m /
        (2.0 * std::numbers::pi * (exit_diameter_m / 2.0));

    return {
        {0.120, 0.280, 0.0},
        {exit_diameter_m, exit_width_m, exit_blockage},
        // The reference MATLAB code uses beta2=-30 deg in its own sign
        // convention. TurboCompGPU uses the complementary positive angle.
        {exit_diameter_m, 60.0, 2.0 * std::sqrt(0.5 * (0.140 * 0.140 + 0.060 * 0.060))}
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 300000.0};
}

OperatingPointRequest request(
    const MapPoint& point,
    const thermomodel::IdealGas& air,
    const ImpellerSystemGeometry& g,
    std::optional<double> slip_override,
    bool rotor_b = false)
{
    // The reference experiment reports total inlet conditions. Convert that
    // state to a static state and inlet velocity on the subsonic ideal-gas
    // branch so the compressor's static inlet boundary has the stated p01/T01.
    constexpr double inlet_total_temperature_k = 288.15;
    constexpr double inlet_total_pressure_pa = 101000.0;
    const double gamma = air.parameters().cp /
        (air.parameters().cp - air.parameters().R);
    const double area1 = effective_flow_area(
        inlet_annulus_area(g.inlet), g.inlet.blockage);
    double mach_lo = 0.0;
    double mach_hi = 1.0;
    for (int iteration = 0; iteration < 80; ++iteration) {
        const double mach = 0.5 * (mach_lo + mach_hi);
        const double temperature = inlet_total_temperature_k /
            (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
        const double pressure = inlet_total_pressure_pa *
            std::pow(temperature / inlet_total_temperature_k,
                     gamma / (gamma - 1.0));
        const auto state = evaluate(air, temperature, pressure);
        const double flow = state.rho * area1 *
            mach * std::sqrt(gamma * air.parameters().R * temperature);
        if (flow < point.mass_flow)
            mach_lo = mach;
        else
            mach_hi = mach;
    }
    const double inlet_mach = 0.5 * (mach_lo + mach_hi);
    const double inlet_temperature_k = inlet_total_temperature_k /
        (1.0 + 0.5 * (gamma - 1.0) * inlet_mach * inlet_mach);
    const double inlet_pressure_pa = inlet_total_pressure_pa *
        std::pow(inlet_temperature_k / inlet_total_temperature_k,
                 gamma / (gamma - 1.0));
    // Map the author's negative-angle convention to TurboCompGPU's
    // positive convention, and use the matching B or A slip-correlation angle.
    const double beta_author_deg = rotor_b ? 40.0 : 30.0;
    constexpr double blades = 20.0;
    const double correlation_slip = 1.0 -
        std::sqrt(std::cos(beta_author_deg * std::numbers::pi / 180.0)) /
        std::pow(blades, 0.7);
    const double slip = slip_override.value_or(correlation_slip);
    if (!(std::isfinite(slip) && slip > 0.0 && slip <= 1.0))
        throw std::invalid_argument("Slip factor override must be in (0,1]");

    return {
        {inlet_temperature_k, inlet_pressure_pa, point.mass_flow,
         point.rpm, slip, 0.0, 0.0},
        {inlet_temperature_k * 1.15,
         std::log(1.5 * inlet_pressure_pa / scales().pressure)}
    };
}

NonlinearSolverOptions options()
{
    NonlinearSolverOptions result;
    result.residual_tolerance = 1.0e-12;
    return result;
}

struct StagePressureRatios {
    double total_to_total;
    double impeller_total_to_total;
    double diffuser_total_pressure_ratio;
};

StagePressureRatios pressure_ratios(
    const thermomodel::IdealGas& air,
    const ImpellerSystemGeometry& g,
    const OperatingPointRequest& op,
    const CompleteStageResult& result)
{
    const auto& station3 = *result.stage->station3;
    const auto& stage = *result.stage;
    const auto state1 = evaluate(
        air, op.operating.inlet_temperature, op.operating.inlet_pressure);
    const double area1 = effective_flow_area(inlet_annulus_area(g.inlet),
                                              g.inlet.blockage);
    const double vm1 = op.operating.mass_flow / (state1.rho * area1);
    const FlowStation station1{state1, {vm1, op.operating.inlet_prewhirl}};
    const FlowStation station2{
        stage.station2_static,
        {stage.station2_meridional_velocity,
         stage.station2_tangential_velocity}};
    const FlowStation station3_flow{
        station3.state,
        {station3.meridional_velocity, station3.tangential_velocity}};
    const double p01 = evaluate_stagnation_state(air, station1)
                           .stagnation_state.p;
    const double p02 = evaluate_stagnation_state(air, station2)
                           .stagnation_state.p;
    const double p03 = evaluate_stagnation_state(air, station3_flow)
                           .stagnation_state.p;
    return {p03 / p01, p02 / p01, p03 / p02};
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path input_path = argc > 1
            ? argv[1]
            : "data/validation/eckardt_impeller_a.csv";
        const std::filesystem::path output_path = argc > 2
            ? argv[2]
            : "benchmarks/results-eckardt-impeller-a/qualification.csv";
        const std::optional<double> slip_override = argc > 3 && argv[3][0] != '\0'
            ? std::optional<double>(std::stod(argv[3]))
            : std::nullopt;
        const bool rotor_b = argc > 4 && std::string(argv[4]) == "B";

        const auto points = read_map(input_path);
        const auto g = geometry(rotor_b);
        const auto s = scales();
        const thermomodel::IdealGas air;

        std::filesystem::create_directories(output_path.parent_path());
        std::ofstream csv(output_path);
        if (!csv)
            throw std::runtime_error("Cannot write: " + output_path.string());
        csv << "rpm,mass_flow_kg_s,slip_factor,"
               "measured_total_pressure_ratio,predicted_total_pressure_ratio,"
               "predicted_impeller_total_pressure_ratio,"
               "diffuser_total_pressure_ratio,relative_error_percent,"
               "stage_complete,message\n";
        csv << std::setprecision(12);

        int complete = 0;
        double sum_error = 0.0;
        double sum_absolute_error = 0.0;
        double max_absolute_error = 0.0;
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();

        for (const auto& point : points) {
            const auto op = request(point, air, g, slip_override, rotor_b);
            const auto result = evaluate_complete_stage(
                air, air, g,
                {g.exit.diameter / 2.0,
                 (rotor_b ? 3.69 : 2.69) * g.exit.diameter / 2.0,
                 0.51 * g.exit.width, 0.0},
                op, s, options());

            double predicted = nan;
            double predicted_impeller = nan;
            double diffuser_ratio = nan;
            double error_percent = nan;
            std::string message = result.message;
            if (result.complete()) {
                const auto ratios = pressure_ratios(air, g, op, result);
                predicted = ratios.total_to_total;
                predicted_impeller = ratios.impeller_total_to_total;
                diffuser_ratio = ratios.diffuser_total_pressure_ratio;
                error_percent = 100.0 *
                    (predicted - point.measured_total_pressure_ratio) /
                    point.measured_total_pressure_ratio;
                ++complete;
                sum_error += error_percent;
                sum_absolute_error += std::abs(error_percent);
                max_absolute_error = std::max(max_absolute_error,
                                              std::abs(error_percent));
            }

            for (char& c : message)
                if (c == ',' || c == '\n' || c == '\r') c = ';';
            csv << point.rpm << ',' << point.mass_flow << ','
                << op.operating.slip_factor << ','
                << point.measured_total_pressure_ratio << ',' << predicted
                << ',' << predicted_impeller << ',' << diffuser_ratio
                << ',' << error_percent << ',' << result.complete() << ','
                << message << '\n';
            std::cout << std::fixed << std::setprecision(4)
                      << point.rpm << " rpm, " << point.mass_flow
                      << " kg/s: measured PR="
                      << point.measured_total_pressure_ratio
                      << ", sigma=" << op.operating.slip_factor;
            if (result.complete())
                std::cout << ", predicted PR=" << predicted
                          << ", error=" << error_percent << "%\n";
            else
                std::cout << ", no complete solution: " << result.message
                          << '\n';
        }

        std::cout << "\nEckardt " << (rotor_b ? "B" : "A")
                  << " points=" << points.size()
                  << " complete=" << complete
                  << " failed=" << points.size() - complete;
        if (complete > 0)
            std::cout << " bias=" << sum_error / complete << "%"
                      << " MAPE=" << sum_absolute_error / complete << "%"
                      << " max_abs=" << max_absolute_error << "%";
        std::cout << "\nCSV: " << output_path << '\n';
        return complete == 0 ? 2 : 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Eckardt validation failed: " << error.what() << '\n';
        return 1;
    }
}
