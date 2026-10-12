#include <turbocompgpu/compressor_map.hpp>

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

constexpr double comparison_tolerance = 1.0e-8;

struct Counts {
    int solves = 0;
    int converged = 0;
    int validated = 0;
    int physical = 0;
    int mapped = 0;
    int paired = 0;
};

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

void check_close(
    const char* name,
    double a,
    double b)
{
    const double scale =
        std::max({1.0, std::abs(a), std::abs(b)});

    if (!std::isfinite(a) ||
        !std::isfinite(b) ||
        std::abs(a - b) > comparison_tolerance * scale) {
        throw std::runtime_error(
            std::string("Paired-map disagreement: ") + name);
    }
}

void compare_maps(
    const CompressorMapPoint& a,
    const CompressorMapPoint& b)
{
    check_close("corrected_mass_flow",
                a.corrected_mass_flow, b.corrected_mass_flow);

    check_close("corrected_speed_rpm",
                a.corrected_speed_rpm, b.corrected_speed_rpm);

    check_close("inlet_total_temperature",
                a.inlet_total_temperature, b.inlet_total_temperature);

    check_close("inlet_total_pressure",
                a.inlet_total_pressure, b.inlet_total_pressure);

    check_close("outlet_total_temperature",
                a.outlet_total_temperature, b.outlet_total_temperature);

    check_close("outlet_total_pressure",
                a.outlet_total_pressure, b.outlet_total_pressure);

    check_close("static_pressure_ratio",
                a.static_pressure_ratio, b.static_pressure_ratio);

    check_close("total_pressure_ratio",
                a.total_pressure_ratio, b.total_pressure_ratio);

    check_close("specific_work",
                a.specific_work, b.specific_work);

    check_close("shaft_power",
                a.shaft_power, b.shaft_power);

    check_close("torque",
                a.torque, b.torque);
}

void write_row(
    std::ofstream& csv,
    const std::string& model,
    const std::string& method,
    const OperatingPointResult& point,
    const CompressorMapPoint& map)
{
    csv << model << ','
        << method << ','
        << point.operating.mass_flow << ','
        << point.operating.rotational_speed_rpm << ','
        << point.converged << ','
        << point.independently_validated << ','
        << point.basic_physics_valid << ','
        << point.iterations << ','
        << point.evaluations << ','
        << point.residual_norm << ','
        << point.validation_residual_norm << ','
        << point.outlet_temperature << ','
        << point.outlet_pressure << ','
        << map.corrected_mass_flow << ','
        << map.corrected_speed_rpm << ','
        << map.inlet_total_temperature << ','
        << map.inlet_total_pressure << ','
        << map.outlet_total_temperature << ','
        << map.outlet_total_pressure << ','
        << map.static_pressure_ratio << ','
        << map.total_pressure_ratio << ','
        << map.specific_work << ','
        << map.shaft_power << ','
        << map.torque << ','
        << map.entropy_generation << ','
        << map.full_residual_norm << ','
        << map.inlet_stagnation_enthalpy_residual << ','
        << map.inlet_stagnation_entropy_residual << ','
        << map.outlet_stagnation_enthalpy_residual << ','
        << map.outlet_stagnation_entropy_residual << ','
        << (map.operating_limits ==
            OperatingLimitStatus::NotEvaluated
                ? "NotEvaluated" : "Evaluated")
        << '\n';
}

CompressorMapPoint evaluate_checked(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& g,
    const NonlinearResidualScales& s,
    const MapReferenceConditions& reference,
    const OperatingPointResult& point,
    Counts& counts)
{
    ++counts.solves;

    if (point.converged)
        ++counts.converged;

    if (point.independently_validated)
        ++counts.validated;

    if (point.basic_physics_valid)
        ++counts.physical;

    if (!point.converged ||
        !point.independently_validated ||
        !point.basic_physics_valid) {
        throw std::runtime_error(
            "Operating-point prerequisite failed");
    }

    const auto map = evaluate_compressor_map_point(
        model, g, point, s, reference);

    if (map.operating_limits !=
        OperatingLimitStatus::NotEvaluated) {
        throw std::runtime_error(
            "Unexpected operating-limit classification");
    }

    check_close(
        "static pressure ratio",
        map.static_pressure_ratio,
        point.pressure_ratio);

    check_close(
        "specific work",
        map.specific_work,
        point.specific_work);

    check_close(
        "shaft power",
        map.shaft_power,
        point.shaft_power);

    ++counts.mapped;
    return map;
}

void sweep(
    const thermomodel::Model& model,
    const std::string& model_name,
    std::ofstream& csv,
    Counts& counts)
{
    const auto g = geometry();
    const auto s = scales();

    const MapReferenceConditions reference {
        288.15, 101325.0
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    constexpr std::array<double, 7> flows {
        2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
    };

    constexpr std::array<double, 5> speeds {
        9000.0, 10500.0, 12000.0, 13500.0, 15000.0
    };

    for (double speed : speeds) {
        NonlinearUnknowns previous {
            340.0, std::log(700000.0 / s.pressure)
        };

        bool previous_valid = false;

        for (double flow : flows) {
            const NonlinearOperatingPoint operating {
                300.0, 500000.0,
                flow, speed,
                0.90, 20.0, 20.0
            };

            const NonlinearUnknowns baseline {
                340.0,
                std::log(700000.0 / s.pressure)
            };

            const auto independent = evaluate_operating_point(
                model, g, {operating, baseline}, s, options);

            const auto continued = evaluate_operating_point(
                model, g,
                {operating, previous_valid ? previous : baseline},
                s, options);

            const auto independent_map = evaluate_checked(
                model, g, s, reference, independent, counts);

            const auto continued_map = evaluate_checked(
                model, g, s, reference, continued, counts);

            compare_maps(independent_map, continued_map);
            ++counts.paired;

            write_row(
                csv, model_name, "independent",
                independent, independent_map);

            write_row(
                csv, model_name, "continuation",
                continued, continued_map);

            previous = {
                continued.outlet_temperature,
                std::log(continued.outlet_pressure / s.pressure)
            };

            previous_valid = true;
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-compressor-map/map.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open CSV output");

        csv << std::setprecision(17);

        csv << "model,method,mass_flow,rotational_speed_rpm,"
               "converged,validated,basic_physics_valid,"
               "iterations,evaluations,residual_norm,"
               "validation_residual_norm,"
               "outlet_temperature,outlet_pressure,"
               "corrected_mass_flow,corrected_speed_rpm,"
               "inlet_total_temperature,inlet_total_pressure,"
               "outlet_total_temperature,outlet_total_pressure,"
               "static_pressure_ratio,total_pressure_ratio,"
               "specific_work,shaft_power,torque,"
               "entropy_generation,full_residual_norm,"
               "inlet_stagnation_enthalpy_residual,"
               "inlet_stagnation_entropy_residual,"
               "outlet_stagnation_enthalpy_residual,"
               "outlet_stagnation_entropy_residual,"
               "operating_limits\n";

        Counts counts;

        const thermomodel::IdealGas ideal;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        const thermomodel::PengRobinson pr(caloric);

        sweep(ideal, "IdealGas", csv, counts);
        sweep(pr, "PengRobinson", csv, counts);

        if (counts.solves != 140 ||
            counts.converged != 140 ||
            counts.validated != 140 ||
            counts.physical != 140 ||
            counts.mapped != 140 ||
            counts.paired != 70) {
            throw std::runtime_error(
                "Compressor-map envelope acceptance criteria failed");
        }

        std::cout
            << "solves=" << counts.solves
            << " converged=" << counts.converged
            << " validated=" << counts.validated
            << " physical=" << counts.physical
            << " mapped=" << counts.mapped
            << " paired=" << counts.paired
            << '\n';

        std::cout << "CSV: " << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
