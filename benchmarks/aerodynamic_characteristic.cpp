#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/incidence_operating_point.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>

#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const NonlinearResidualScales scales{
    50000.0, 100.0, 500000.0
};

const MapReferenceConditions reference{
    288.15, 101325.0
};

struct Counts {
    int attempted = 0;
    int mapped = 0;
    int failed = 0;
};

void write_row(
    std::ofstream& csv,
    const std::string& model_name,
    double speed,
    double flow,
    bool success,
    const OperatingPointResult* point,
    const CompressorMapPoint* map,
    double incidence_entropy,
    double diffusion_entropy,
    double incidence_loss,
    double diffusion_loss)
{
    const double nan =
        std::numeric_limits<double>::quiet_NaN();

    csv
        << model_name << ','
        << speed << ','
        << flow << ','
        << success << ','
        << (point ? point->converged : false) << ','
        << (point ? point->independently_validated : false) << ','
        << (point ? point->basic_physics_valid : false) << ','
        << (map ? map->total_pressure_ratio : nan) << ','
        << (map ? map->static_pressure_ratio : nan) << ','
        << (map ? map->shaft_power : nan) << ','
        << (map ? map->corrected_mass_flow : nan) << ','
        << (map ? map->corrected_speed_rpm : nan) << ','
        << (point ? point->outlet_temperature : nan) << ','
        << (point ? point->outlet_pressure : nan) << ','
        << (point ? point->entropy_generation : nan) << ','
        << incidence_entropy << ','
        << diffusion_entropy << ','
        << incidence_loss << ','
        << diffusion_loss << ','
        << (point ? point->validation_residual_norm : nan)
        << '\n';
}

void require_valid(
    const OperatingPointResult& point,
    const std::string& label)
{
    if (!(point.converged &&
          point.independently_validated &&
          point.basic_physics_valid)) {
        throw std::runtime_error(
            label + ": " + point.message);
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-aerodynamic-characteristic/"
                  "comparison.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open CSV");

        csv << std::setprecision(17);

        csv
            << "model,rotational_speed_rpm,mass_flow,"
            << "mapped,converged,validated,physical,"
            << "total_pressure_ratio,static_pressure_ratio,"
            << "shaft_power,corrected_mass_flow,"
            << "corrected_speed_rpm,"
            << "outlet_temperature,outlet_pressure,"
            << "total_entropy_generation,"
            << "incidence_entropy_generation,"
            << "diffusion_entropy_generation,"
            << "incidence_specific_loss,"
            << "diffusion_specific_loss,"
            << "validation_residual_norm\n";

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane(
            caloric, parameters);

        const auto inlet = methane.evaluate_Tp(
            300.0, 500000.0);

        if (!inlet)
            throw std::runtime_error(
                "Inlet EOS evaluation failed");

        // Preserve the established reference calibration.
        const OperatingIncidenceConditions design{
            5.0,
            12000.0,
            inlet.value.rho,
            20.0
        };

        const double blade_angle =
            calibrate_inlet_blade_angle_deg(
                geometry, design);

        const IncidenceLossParameters incidence_params{
            blade_angle, 1.0
        };

        const DiffusionLossParameters diffusion_params{
            1.0, 1.5
        };

        NonlinearSolverOptions options;
        options.residual_tolerance = 1.0e-13;

        Counts counts;

        std::cout << std::setprecision(12)
                  << "blade_angle_deg=" << blade_angle
                  << '\n';

        for (double speed : speeds) {
            for (int i = 0; i <= 20; ++i) {
                const double flow = 0.5 + 0.1 * i;

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

                const OperatingPointRequest request{
                    operating, guess
                };

                for (int variant = 0; variant < 3; ++variant) {
                    const std::string name =
                        variant == 0 ? "baseline" :
                        variant == 1 ? "incidence" :
                                       "combined";

                    ++counts.attempted;

                    try {
                        OperatingPointResult point;

                        double incidence_entropy = 0.0;
                        double diffusion_entropy = 0.0;
                        double incidence_loss = 0.0;
                        double diffusion_loss = 0.0;

                        if (variant == 0) {
                            point = evaluate_operating_point(
                                methane,
                                geometry,
                                request,
                                scales,
                                options);
                        }
                        else if (variant == 1) {
                            const auto result =
                                evaluate_incidence_operating_point(
                                    methane,
                                    geometry,
                                    request,
                                    scales,
                                    incidence_params,
                                    options);

                            point = result.solution;

                            incidence_entropy =
                                result.incidence_entropy_generation;

                            incidence_loss =
                                result.incidence.loss.specific_loss;
                        }
                        else {
                            const auto result =
                                evaluate_diffusion_operating_point(
                                    methane,
                                    geometry,
                                    request,
                                    scales,
                                    incidence_params,
                                    diffusion_params,
                                    options);

                            point = result.validation;

                            if (!result.independently_validated)
                                throw std::runtime_error(
                                    "Coupled diffusion validation failed");

                            incidence_entropy =
                                result.incidence_entropy_generation;

                            diffusion_entropy =
                                result.diffusion_entropy_generation;

                            incidence_loss =
                                result.incidence.loss.specific_loss;

                            diffusion_loss =
                                result.diffusion.loss.specific_loss;
                        }

                        require_valid(point, name);

                        const auto map =
                            evaluate_compressor_map_point(
                                methane,
                                geometry,
                                point,
                                scales,
                                reference);

                        write_row(
                            csv, name, speed, flow, true,
                            &point, &map,
                            incidence_entropy,
                            diffusion_entropy,
                            incidence_loss,
                            diffusion_loss);

                        ++counts.mapped;
                    }
                    catch (const std::exception& error) {
                        write_row(
                            csv, name, speed, flow, false,
                            nullptr, nullptr,
                            0.0, 0.0, 0.0, 0.0);

                        ++counts.failed;

                        std::cerr
                            << "FAILED"
                            << " model=" << name
                            << " rpm=" << speed
                            << " flow=" << flow
                            << " reason=" << error.what()
                            << '\n';
                    }
                }
            }
        }

        std::cout
            << "\nSTAGE 7E SUMMARY\n"
            << "attempted=" << counts.attempted << '\n'
            << "mapped=" << counts.mapped << '\n'
            << "failed=" << counts.failed << '\n'
            << "output=" << output << '\n';

        if (counts.attempted != 315)
            throw std::runtime_error(
                "Incorrect number of attempted cases");

        return counts.failed == 0 ? 0 : 2;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
