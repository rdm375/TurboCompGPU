#include <turbocompgpu/transient_pipe_coupling.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <pipe_sim/physics.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

struct Config {
    int cells = 12;
    double pipe_length_m = 1200.0;
    double dt_s = 1.0;
    double duration_s = 30.0;
    double theta = 0.65;
    double diameter_m = 0.60;
    double temperature_k = 300.0;
    double compressibility = 1.0;
    double gas_constant = 287.05;
    double viscosity_pa_s = 1.8e-5;
    double roughness_m = 1.0e-5;
    double source_pressure_pa = 500000.0;
    double initial_flow_kg_s = 5.0;
    double export_flow_before_kg_s = 5.0;
    double export_flow_after_kg_s = 5.5;
    double demand_change_time_s = 10.0;
    double speed_before_rpm = 12000.0;
    double speed_after_rpm = 12000.0;
    double speed_change_time_s = 15.0;
    double minimum_throughflow_kg_s = 3.0;
    double maximum_throughflow_kg_s = 8.0;
    double slip_factor = 0.90;
    double inlet_hub_diameter_m = 0.10;
    double inlet_shroud_diameter_m = 0.20;
    double inlet_blockage = 0.05;
    double impeller_inlet_diameter_m = 0.15;
    double impeller_exit_diameter_m = 0.40;
    double impeller_exit_width_m = 0.025;
    double impeller_exit_blockage = 0.10;
    double blade_angle_deg = 45.0;
    double impeller_inlet_blockage = 0.15;
    double diffuser_radius_ratio = 1.50;
};

double number(const std::string& value, const std::string& key)
{
    std::size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result))
        throw std::runtime_error("Invalid numeric value for " + key);
    return result;
}

Config read_config(const std::filesystem::path& path)
{
    Config c;
    if (path.empty()) return c;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open config: " + path.string());
    std::map<std::string, double*> fields{
        {"pipe_length_m", &c.pipe_length_m}, {"dt_s", &c.dt_s},
        {"duration_s", &c.duration_s}, {"theta", &c.theta},
        {"diameter_m", &c.diameter_m}, {"temperature_k", &c.temperature_k},
        {"compressibility", &c.compressibility}, {"gas_constant", &c.gas_constant},
        {"viscosity_pa_s", &c.viscosity_pa_s}, {"roughness_m", &c.roughness_m},
        {"source_pressure_pa", &c.source_pressure_pa},
        {"initial_flow_kg_s", &c.initial_flow_kg_s},
        {"export_flow_before_kg_s", &c.export_flow_before_kg_s},
        {"export_flow_after_kg_s", &c.export_flow_after_kg_s},
        {"demand_change_time_s", &c.demand_change_time_s},
        {"speed_before_rpm", &c.speed_before_rpm},
        {"speed_after_rpm", &c.speed_after_rpm},
        {"speed_change_time_s", &c.speed_change_time_s},
        {"minimum_throughflow_kg_s", &c.minimum_throughflow_kg_s},
        {"maximum_throughflow_kg_s", &c.maximum_throughflow_kg_s},
        {"slip_factor", &c.slip_factor},
        {"inlet_hub_diameter_m", &c.inlet_hub_diameter_m},
        {"inlet_shroud_diameter_m", &c.inlet_shroud_diameter_m},
        {"inlet_blockage", &c.inlet_blockage},
        {"impeller_inlet_diameter_m", &c.impeller_inlet_diameter_m},
        {"impeller_exit_diameter_m", &c.impeller_exit_diameter_m},
        {"impeller_exit_width_m", &c.impeller_exit_width_m},
        {"impeller_exit_blockage", &c.impeller_exit_blockage},
        {"blade_angle_deg", &c.blade_angle_deg},
        {"impeller_inlet_blockage", &c.impeller_inlet_blockage},
        {"diffuser_radius_ratio", &c.diffuser_radius_ratio},
    };
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.resize(comment);
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        const auto equal = line.find('=', first);
        if (equal == std::string::npos)
            throw std::runtime_error("Expected key=value at config line " +
                                     std::to_string(line_number));
        const auto last = line.find_last_not_of(" \t\r", equal - 1);
        const std::string key = line.substr(first, last - first + 1);
        const auto value_start = line.find_first_not_of(" \t", equal + 1);
        if (value_start == std::string::npos)
            throw std::runtime_error("Missing value for " + key);
        const auto value_end = line.find_last_not_of(" \t\r");
        if (key == "cells") {
            const double parsed = number(line.substr(value_start, value_end - value_start + 1), key);
            if (parsed < 2.0 || parsed > 10000.0 || std::floor(parsed) != parsed)
                throw std::runtime_error("cells must be an integer in [2, 10000]");
            c.cells = static_cast<int>(parsed);
            continue;
        }
        const auto field = fields.find(key);
        if (field == fields.end()) throw std::runtime_error("Unknown config key: " + key);
        *field->second = number(line.substr(value_start, value_end - value_start + 1), key);
    }
    if (c.cells < 2 || !(c.pipe_length_m > 0.0 && c.dt_s > 0.0 &&
        c.duration_s > 0.0 && c.theta >= 0.0 && c.theta <= 1.0 &&
        c.diameter_m > 0.0 && c.temperature_k > 0.0 &&
        c.compressibility > 0.0 && c.gas_constant > 0.0 &&
        c.viscosity_pa_s > 0.0 && c.roughness_m >= 0.0 &&
        c.source_pressure_pa > 0.0 && c.initial_flow_kg_s > 0.0 &&
        c.minimum_throughflow_kg_s > 0.0 &&
        c.maximum_throughflow_kg_s > c.minimum_throughflow_kg_s &&
        c.speed_before_rpm > 0.0 && c.speed_after_rpm > 0.0 &&
        c.slip_factor > 0.0 && c.slip_factor <= 1.0 &&
        c.inlet_hub_diameter_m >= 0.0 &&
        c.inlet_shroud_diameter_m > c.inlet_hub_diameter_m &&
        c.inlet_blockage >= 0.0 && c.inlet_blockage < 1.0 &&
        c.impeller_inlet_diameter_m > 0.0 && c.impeller_exit_diameter_m > 0.0 &&
        c.impeller_exit_width_m > 0.0 && c.impeller_exit_blockage >= 0.0 &&
        c.impeller_exit_blockage < 1.0 && c.diffuser_radius_ratio > 1.0))
        throw std::runtime_error("Config has an out-of-range parameter");
    if (!(c.impeller_inlet_blockage >= 0.0 && c.impeller_inlet_blockage < 1.0))
        throw std::runtime_error("impeller_inlet_blockage must be in [0,1)");
    return c;
}

PipeSegmentModel make_pipe(const Config& c)
{
    const double area = std::numbers::pi * c.diameter_m * c.diameter_m / 4.0;
    return {c.cells,
        {c.pipe_length_m / c.cells, c.dt_s, c.theta, c.diameter_m, area,
         c.temperature_k, c.compressibility, c.gas_constant,
         c.viscosity_pa_s, c.roughness_m},
        {1.0e-8, 1.0e-12, 60, false}};
}

PipeSegmentState initialize_pipe(
    const PipeSegmentModel& model, double inlet_pressure_pa, double flow_kg_s)
{
    const auto& p = model.parameters;
    const double reynolds = pipe_sim::reynolds_mass(
        flow_kg_s, p.diameter, p.area, p.viscosity);
    const double darcy = pipe_sim::friction_sj(
        reynolds, p.roughness, p.diameter);
    const double coefficient = darcy * p.z * p.gas_constant * p.temperature *
        flow_kg_s * flow_kg_s / (p.diameter * p.area * p.area);
    PipeSegmentState state;
    state.pressure_pa.resize(static_cast<std::size_t>(model.cells + 1));
    state.flow_kg_s.assign(static_cast<std::size_t>(model.cells), flow_kg_s);
    state.inlet_flow_kg_s = state.outlet_flow_kg_s = flow_kg_s;
    for (int i = 0; i <= model.cells; ++i) {
        const double p2 = inlet_pressure_pa * inlet_pressure_pa -
            coefficient * static_cast<double>(i) * p.dx;
        if (!(p2 > 0.0)) throw std::runtime_error("Pipe initializer has nonpositive pressure squared");
        state.pressure_pa[static_cast<std::size_t>(i)] = std::sqrt(p2);
    }
    return state;
}

void write_leg(
    std::ofstream& out, double time_s, const char* leg, double offset_m,
    const PipeSegmentModel& model, const PipeSegmentState& state)
{
    const int n = model.cells;
    for (int i = 0; i <= n; ++i) {
        double flow;
        if (i == 0) flow = state.inlet_flow_kg_s;
        else if (i == n) flow = state.outlet_flow_kg_s;
        else flow = 0.5 * (state.flow_kg_s[static_cast<std::size_t>(i - 1)] +
                            state.flow_kg_s[static_cast<std::size_t>(i)]);
        // A one-metre plotting gap marks the compressor station between legs.
        const double x = offset_m + model.parameters.dx * i;
        out << time_s << ',' << leg << ',' << x << ','
            << state.pressure_pa[static_cast<std::size_t>(i)] << ','
            << flow << ',' << model.parameters.temperature << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::filesystem::path config_path;
        std::filesystem::path output_dir = "two_pipe_compressor_output";
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--config" && i + 1 < argc) config_path = argv[++i];
            else if (arg == "--output-dir" && i + 1 < argc) output_dir = argv[++i];
            else throw std::runtime_error(
                "Usage: turbocompgpu_two_pipe_compressor [--config FILE] [--output-dir DIR]");
        }
        const Config c = read_config(config_path);
        const double step_count = c.duration_s / c.dt_s;
        if (std::abs(step_count - std::round(step_count)) > 1.0e-10)
            throw std::runtime_error("duration_s must be an integer multiple of dt_s");
        std::filesystem::create_directories(output_dir);

        thermomodel::IdealGas gas;
        ImpellerSystemGeometry impeller{
            {c.inlet_hub_diameter_m, c.inlet_shroud_diameter_m, c.inlet_blockage},
            {c.impeller_exit_diameter_m, c.impeller_exit_width_m,
             c.impeller_exit_blockage},
            {c.impeller_exit_diameter_m, c.blade_angle_deg,
             c.impeller_inlet_diameter_m}};
        const double r2 = impeller.exit.diameter / 2.0;
        const double exit_area = effective_flow_area(
            impeller_exit_area(impeller.exit), impeller.exit.blockage);
        const VanelessDiffuserGeometry diffuser{
            r2, c.diffuser_radius_ratio * r2,
            exit_area / (2.0 * std::numbers::pi * r2), 0.0};
        CompressorComponentConfig compressor_config;
        compressor_config.slip_factor = c.slip_factor;
        compressor_config.residual_scales = {50000.0, 100.0, 500000.0};
        compressor_config.solver_options.residual_tolerance = 1.0e-13;
        compressor_config.initial_outlet_temperature_k = c.temperature_k + 40.0;
        compressor_config.initial_outlet_pressure_pa = 1.4 * c.source_pressure_pa;
        const CompressorComponent compressor(
            gas, gas, impeller, diffuser, compressor_config);
        const auto pipe = make_pipe(c);
        auto upstream = initialize_pipe(pipe, c.source_pressure_pa, c.initial_flow_kg_s);
        const auto initial_compressor = compressor.evaluate({
            {c.temperature_k, upstream.pressure_pa.back()},
            c.initial_flow_kg_s, c.speed_before_rpm, std::nullopt});
        if (!initial_compressor.valid || !initial_compressor.outlet_static_state)
            throw std::runtime_error("Initial compressor point failed: " + initial_compressor.message);
        auto downstream = initialize_pipe(
            pipe, initial_compressor.outlet_static_state->p, c.initial_flow_kg_s);

        std::ofstream field(output_dir / "pipeline_fields.csv");
        std::ofstream operation(output_dir / "compressor_history.csv");
        if (!field || !operation) throw std::runtime_error("Cannot create output CSV files");
        field << std::setprecision(12)
              << "time_s,leg,distance_m,pressure_pa,flow_kg_s,temperature_k\n";
        operation << "time_s,source_pressure_pa,suction_pressure_pa,discharge_pressure_pa,"
                     "export_flow_kg_s,compressor_throughflow_kg_s,shaft_power_w,"
                     "compressor_node_mass_residual_kg_s,speed_rpm\n";
        const double length_upstream = pipe.parameters.dx * pipe.cells;
        write_leg(field, 0.0, "upstream", 0.0, pipe, upstream);
        write_leg(field, 0.0, "downstream", length_upstream + 1.0, pipe, downstream);
        TwoPipeCompressorStepOptions options;
        options.minimum_compressor_throughflow_kg_s = c.minimum_throughflow_kg_s;
        options.maximum_compressor_throughflow_kg_s = c.maximum_throughflow_kg_s;
        const int steps = static_cast<int>(std::ceil(c.duration_s / c.dt_s));
        for (int k = 1; k <= steps; ++k) {
            const double time = std::min(k * c.dt_s, c.duration_s);
            const double export_flow = time < c.demand_change_time_s
                ? c.export_flow_before_kg_s : c.export_flow_after_kg_s;
            const double speed = time < c.speed_change_time_s
                ? c.speed_before_rpm : c.speed_after_rpm;
            const auto result = advance_two_pipe_compressor(
                pipe, upstream, c.source_pressure_pa, pipe, downstream,
                export_flow, compressor, speed, options);
            if (!result.valid) throw std::runtime_error(
                "Coupled step failed at t=" + std::to_string(time) + " s: " + result.message);
            operation << std::setprecision(12) << time << ',' << c.source_pressure_pa << ','
                      << result.measurements.suction_pressure_pa << ','
                      << result.measurements.discharge_pressure_pa << ','
                      << result.measurements.pipeline_export_flow_kg_s << ','
                      << result.compressor.compressor_throughflow_kg_s << ','
                      << result.measurements.compressor_shaft_power_w << ','
                      << result.compressor_node_mass_residual_kg_s << ',' << speed << '\n';
            write_leg(field, time, "upstream", 0.0, pipe, result.upstream);
            write_leg(field, time, "downstream", length_upstream + 1.0, pipe, result.downstream);
            upstream = result.upstream;
            downstream = result.downstream;
        }
        std::cout << "Wrote " << (output_dir / "pipeline_fields.csv") << "\n"
                  << "Wrote " << (output_dir / "compressor_history.csv") << "\n"
                  << "Isothermal pipe temperature is fixed per segment; temperature maps show those configured values.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "two_pipe_compressor: " << error.what() << '\n';
        return 1;
    }
}
