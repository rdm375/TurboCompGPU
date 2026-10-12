#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/passage_turning_choking.hpp>

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

using namespace turbocompgpu;

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-choking-characteristic/"
                  "endpoints.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open output CSV");

        csv << std::setprecision(17);

        csv
            << "rotational_speed_rpm,"
            << "critical_mass_flow,"
            << "total_pressure_ratio,"
            << "static_pressure_ratio,"
            << "shaft_power,"
            << "corrected_mass_flow,"
            << "corrected_speed_rpm,"
            << "inlet_static_temperature,"
            << "inlet_static_pressure,"
            << "inlet_total_temperature,"
            << "inlet_total_pressure,"
            << "inlet_total_temperature_difference,"
            << "inlet_total_pressure_difference,"
            << "outlet_total_temperature,"
            << "outlet_total_pressure,"
            << "converged,"
            << "validated,"
            << "basic_physics_valid,"
            << "operating_limits_evaluated,"
            << "residual_norm,"
            << "validation_residual_norm,"
            << "choking_mass_residual,"
            << "choking_rothalpy_residual,"
            << "choking_continuity_residual,"
            << "choking_sonic_residual\n";

        const ImpellerSystemGeometry geometry{
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {0.40, 45.0, 0.15}
        };

        const ThroatGeometry throat{
            12, 0.15, 0.0035, 45.0, 0.002, 0.10
        };

        const NonlinearResidualScales scales{
            50000.0, 100.0, 500000.0
        };

        const MapReferenceConditions reference{
            288.15, 101325.0
        };

        constexpr std::array<double, 5> speeds{
            9000.0,
            10500.0,
            12000.0,
            13500.0,
            15000.0
        };

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        thermomodel::PengRobinson methane(
            caloric, parameters);

        NonlinearSolverOptions options;
        options.residual_tolerance = 1.0e-13;

        int validated_endpoints = 0;

        for (double speed : speeds) {
            const PassageTurningChokingInputs choking_input{
                300.0,
                500000.0,
                speed,
                20.0,
                throat,
                0.5
            };

            const auto choking =
                evaluate_passage_turning_choking(
                    methane,
                    methane,
                    geometry,
                    choking_input);

            const double critical_flow =
                choking.choking.coupled.critical_mass_flow;

            if (!(std::isfinite(critical_flow) &&
                  critical_flow > 0.0))
                throw std::runtime_error(
                    "Invalid critical mass flow");

            const NonlinearOperatingPoint operating{
                300.0,
                500000.0,
                critical_flow,
                speed,
                0.90,
                20.0,
                20.0
            };

            const NonlinearUnknowns guess{
                340.0,
                std::log(700000.0 / scales.pressure)
            };

            const auto point =
                evaluate_operating_point(
                    methane,
                    geometry,
                    {operating, guess},
                    scales,
                    options);

            if (!(point.converged &&
                  point.independently_validated &&
                  point.basic_physics_valid)) {
                throw std::runtime_error(
                    "Endpoint operating-point validation failed "
                    "at RPM=" + std::to_string(speed) +
                    ": " + point.message);
            }

            const auto map =
                evaluate_compressor_map_point(
                    methane,
                    geometry,
                    point,
                    scales,
                    reference);

            // The operating-point solver does not yet evaluate
            // choking or surge limits.
            if (map.operating_limits !=
                OperatingLimitStatus::NotEvaluated) {
                throw std::runtime_error(
                    "Unexpected operating-limit classification");
            }

            ++validated_endpoints;

            const double delta_T0 =
                map.inlet_total_temperature -
                operating.inlet_temperature;

            const double delta_p0 =
                map.inlet_total_pressure -
                operating.inlet_pressure;

            csv
                << speed << ','
                << critical_flow << ','
                << map.total_pressure_ratio << ','
                << map.static_pressure_ratio << ','
                << map.shaft_power << ','
                << map.corrected_mass_flow << ','
                << map.corrected_speed_rpm << ','
                << operating.inlet_temperature << ','
                << operating.inlet_pressure << ','
                << map.inlet_total_temperature << ','
                << map.inlet_total_pressure << ','
                << delta_T0 << ','
                << delta_p0 << ','
                << map.outlet_total_temperature << ','
                << map.outlet_total_pressure << ','
                << point.converged << ','
                << point.independently_validated << ','
                << point.basic_physics_valid << ','
                << 0 << ','
                << point.residual_norm << ','
                << point.validation_residual_norm << ','
                << choking.mass_flow_residual << ','
                << choking.choking.rothalpy_residual << ','
                << choking.choking.continuity_residual << ','
                << choking.choking.sonic_residual
                << '\n';

            std::cout
                << std::setprecision(12)
                << "RPM=" << speed
                << " mdot_critical=" << critical_flow
                << " Pi_tt=" << map.total_pressure_ratio
                << " power_kW=" << map.shaft_power / 1000.0
                << " delta_T0=" << delta_T0
                << " delta_p0=" << delta_p0
                << " validated=1"
                << '\n';
        }

        if (validated_endpoints != 5)
            throw std::runtime_error(
                "Incorrect number of validated endpoints");

        std::cout
            << "\nSTAGE 7D SUMMARY\n"
            << "endpoints=" << validated_endpoints << '\n'
            << "all_validated=1\n"
            << "choke_admissibility_established=0\n"
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
