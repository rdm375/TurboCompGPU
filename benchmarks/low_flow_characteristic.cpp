#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/passage_turning_choking.hpp>

#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

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

constexpr double nan =
    std::numeric_limits<double>::quiet_NaN();

struct Counts {
    int attempted = 0;
    int converged = 0;
    int validated = 0;
    int physical = 0;
    int mapped = 0;
    int failed = 0;
    int above_capacity = 0;
};

void write_row(
    std::ofstream& csv,
    double speed,
    double flow,
    double critical_flow,
    const OperatingPointResult& point,
    const CompressorMapPoint* map)
{
    const double margin =
        1.0 - flow / critical_flow;

    const bool above_capacity = margin < 0.0;

    csv
        << speed << ','
        << flow << ','
        << critical_flow << ','
        << margin << ','
        << static_cast<int>(above_capacity) << ','
        << static_cast<int>(point.converged) << ','
        << static_cast<int>(point.independently_validated) << ','
        << static_cast<int>(point.basic_physics_valid) << ','
        << static_cast<int>(map != nullptr) << ','
        << static_cast<int>(point.status) << ','
        << point.iterations << ','
        << point.evaluations << ','
        << point.residual_norm << ','
        << point.validation_residual_norm << ','
        << point.outlet_temperature << ','
        << point.outlet_pressure << ','
        << point.specific_work << ','
        << point.shaft_power << ','
        << point.torque << ','
        << point.entropy_generation << ',';

    if (map) {
        csv
            << map->corrected_mass_flow << ','
            << map->corrected_speed_rpm << ','
            << map->static_pressure_ratio << ','
            << map->total_pressure_ratio << ','
            << map->inlet_total_temperature << ','
            << map->inlet_total_pressure << ','
            << map->outlet_total_temperature << ','
            << map->outlet_total_pressure;
    } else {
        csv
            << nan << ','
            << nan << ','
            << nan << ','
            << nan << ','
            << nan << ','
            << nan << ','
            << nan << ','
            << nan;
    }

    csv << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-low-flow-characteristic/"
                  "characteristic.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open CSV");

        csv << std::setprecision(17);

        csv
            << "rotational_speed_rpm,mass_flow,"
            << "critical_mass_flow,capacity_margin,"
            << "above_independent_capacity,"
            << "converged,validated,basic_physics_valid,"
            << "mapped,status,iterations,evaluations,"
            << "residual_norm,validation_residual_norm,"
            << "outlet_temperature,outlet_pressure,"
            << "specific_work,shaft_power,torque,"
            << "entropy_generation,"
            << "corrected_mass_flow,corrected_speed_rpm,"
            << "static_pressure_ratio,total_pressure_ratio,"
            << "inlet_total_temperature,inlet_total_pressure,"
            << "outlet_total_temperature,outlet_total_pressure\n";

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

        Counts counts;

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
                  critical_flow > 0.0)) {
                throw std::runtime_error(
                    "Invalid choking capacity");
            }

            std::cout
                << "RPM=" << speed
                << " critical_flow=" << critical_flow
                << '\n';

            for (int i = 0; i <= 20; ++i) {
                const double flow =
                    0.5 + 0.1 * static_cast<double>(i);

                const NonlinearOperatingPoint operating{
                    300.0,
                    500000.0,
                    flow,
                    speed,
                    0.90,
                    20.0,
                    20.0
                };

                const NonlinearUnknowns guess{
                    340.0,
                    std::log(700000.0 / scales.pressure)
                };

                ++counts.attempted;

                OperatingPointResult point;
                point.operating = operating;

                try {
                    point = evaluate_operating_point(
                        methane,
                        geometry,
                        {operating, guess},
                        scales,
                        options);
                }
                catch (const std::exception& error) {
                    std::cerr
                        << "Operating-point exception: "
                        << error.what() << '\n';
                }

                if (point.converged)
                    ++counts.converged;

                if (point.independently_validated)
                    ++counts.validated;

                if (point.basic_physics_valid)
                    ++counts.physical;

                if (flow > critical_flow)
                    ++counts.above_capacity;

                bool mapped = false;

                if (point.converged &&
                    point.independently_validated &&
                    point.basic_physics_valid) {

                    try {
                        const auto map =
                            evaluate_compressor_map_point(
                                methane,
                                geometry,
                                point,
                                scales,
                                reference);

                        write_row(
                            csv,
                            speed,
                            flow,
                            critical_flow,
                            point,
                            &map);

                        ++counts.mapped;
                        mapped = true;
                    }
                    catch (const std::exception& error) {
                        std::cerr
                            << "Map evaluation failed at "
                            << speed << " RPM, "
                            << flow << " kg/s: "
                            << error.what() << '\n';
                    }
                }

                if (!mapped) {
                    ++counts.failed;

                    write_row(
                        csv,
                        speed,
                        flow,
                        critical_flow,
                        point,
                        nullptr);
                }

                std::cout
                    << "  flow=" << flow
                    << " converged=" << point.converged
                    << " validated="
                    << point.independently_validated
                    << " physical="
                    << point.basic_physics_valid
                    << " mapped=" << mapped
                    << " margin="
                    << 1.0 - flow / critical_flow
                    << '\n';
            }
        }

        if (counts.attempted != 105 ||
            counts.mapped + counts.failed != 105) {
            throw std::runtime_error(
                "Low-flow sweep accounting failed");
        }

        std::cout
            << "\nSTAGE 7B SUMMARY\n"
            << "attempted=" << counts.attempted << '\n'
            << "converged=" << counts.converged << '\n'
            << "validated=" << counts.validated << '\n'
            << "physical=" << counts.physical << '\n'
            << "mapped=" << counts.mapped << '\n'
            << "failed=" << counts.failed << '\n'
            << "above_independent_capacity="
            << counts.above_capacity << '\n'
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
