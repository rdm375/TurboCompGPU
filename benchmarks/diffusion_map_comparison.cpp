#include <turbocompgpu/incidence_operating_point.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>
#include <turbocompgpu/compressor_map.hpp>

#include <thermomodel/ideal_gas.hpp>
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

namespace {

using namespace turbocompgpu;

constexpr double Ru = 8.31446261815324;
constexpr double methane_cp = 2200.0;
constexpr double incidence_coefficient = 1.0;

const ImpellerSystemGeometry geometry {
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const NonlinearResidualScales scales {
    50000.0, 100.0, 500000.0
};

const MapReferenceConditions reference {
    288.15, 101325.0
};

constexpr std::array<double, 7> flows {
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
};

constexpr std::array<double, 5> speeds {
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

struct Counts {
    int points = 0;
    int baseline = 0;
    int zero = 0;
    int incidence = 0;
    int mapped = 0;
    int combined = 0;
    int zero_diffusion = 0;
    int active_diffusion = 0;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

bool close(double a, double b, double tolerance = 1e-10)
{
    return std::isfinite(a) &&
           std::isfinite(b) &&
           std::abs(a - b) <= tolerance *
               std::max({1.0, std::abs(a), std::abs(b)});
}

void require_valid(
    const OperatingPointResult& point,
    const std::string& label)
{
    require(
        point.converged &&
        point.independently_validated &&
        point.basic_physics_valid,
        label + ": " + point.message
    );
}

void require_equal(
    const OperatingPointResult& baseline,
    const OperatingPointResult& zero)
{
    require(
        close(baseline.outlet_temperature,
              zero.outlet_temperature) &&
        close(baseline.outlet_pressure,
              zero.outlet_pressure) &&
        close(baseline.specific_work,
              zero.specific_work) &&
        close(baseline.shaft_power,
              zero.shaft_power) &&
        close(baseline.torque,
              zero.torque),
        "Zero-coefficient regression mismatch"
    );
}

void write_header(std::ofstream& csv)
{
    csv
        << "model,mass_flow,rotational_speed_rpm,"
        << "blade_angle_deg,incidence_angle_deg,"
        << "inlet_density,inlet_relative_velocity,"
        << "incidence_specific_loss,"
        << "baseline_entropy_generation,"
        << "incidence_entropy_generation,"
        << "total_entropy_generation,"
        << "baseline_outlet_temperature,"
        << "incidence_outlet_temperature,"
        << "delta_outlet_temperature,"
        << "baseline_outlet_pressure,"
        << "incidence_outlet_pressure,"
        << "delta_outlet_pressure,"
        << "baseline_static_pressure_ratio,"
        << "incidence_static_pressure_ratio,"
        << "delta_static_pressure_ratio,"
        << "baseline_total_pressure_ratio,"
        << "incidence_total_pressure_ratio,"
        << "delta_total_pressure_ratio,"
        << "baseline_specific_work,"
        << "incidence_specific_work,"
        << "delta_specific_work,"
        << "baseline_shaft_power,"
        << "incidence_shaft_power,"
        << "delta_shaft_power,"
        << "baseline_torque,"
        << "incidence_torque,"
        << "delta_torque,"
        << "baseline_corrected_mass_flow,"
        << "incidence_corrected_mass_flow,"
        << "baseline_corrected_speed_rpm,"
        << "incidence_corrected_speed_rpm,"
        << "baseline_full_residual_norm,"
        << "incidence_full_residual_norm,"
        << "diffusion_ratio,"
        << "diffusion_specific_loss,"
        << "diffusion_entropy_generation,"
        << "combined_total_entropy_generation,"
        << "combined_static_pressure_ratio,"
        << "combined_total_pressure_ratio,"
        << "combined_outlet_temperature,"
        << "combined_specific_work,"
        << "combined_shaft_power,"
        << "combined_full_residual_norm\n";
}

void sweep(
    const thermomodel::Model& model,
    const std::string& name,
    std::ofstream& csv,
    Counts& counts)
{
    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto inlet = model.evaluate_Tp(
        300.0, 500000.0
    );

    require(
        static_cast<bool>(inlet),
        name + ": inlet thermodynamic evaluation failed"
    );

    // Calibrate once per thermodynamic model.
    // The blade angle is fixed throughout the sweep.
    const OperatingIncidenceConditions design {
        5.0,
        12000.0,
        inlet.value.rho,
        20.0
    };

    const double blade_angle =
        calibrate_inlet_blade_angle_deg(
            geometry, design
        );

    const IncidenceLossParameters active_params {
        blade_angle,
        incidence_coefficient
    };

    const IncidenceLossParameters zero_params {
        blade_angle,
        0.0
    };

    const DiffusionLossParameters diffusion_params {
        1.0, 1.5
    };

    const DiffusionLossParameters zero_diffusion_params {
        0.0, 1.5
    };

    std::cout
        << name
        << " blade_angle_deg=" << blade_angle
        << '\n';

    for (double speed : speeds) {
        for (double flow : flows) {
            const NonlinearOperatingPoint operating {
                300.0,
                500000.0,
                flow,
                speed,
                0.90,
                20.0,
                20.0
            };

            const NonlinearUnknowns guess {
                340.0,
                std::log(700000.0 / scales.pressure)
            };

            const OperatingPointRequest request {
                operating,
                guess
            };

            const std::string label =
                name + " flow=" + std::to_string(flow) +
                " speed=" + std::to_string(speed);

            ++counts.points;

            const auto baseline =
                evaluate_operating_point(
                    model, geometry, request, scales, options
                );

            require_valid(baseline, label + " baseline");
            ++counts.baseline;

            const auto zero =
                evaluate_incidence_operating_point(
                    model, geometry, request, scales,
                    zero_params, options
                );

            require_valid(
                zero.solution, label + " zero coefficient"
            );

            require_equal(baseline, zero.solution);
            ++counts.zero;

            const auto active =
                evaluate_incidence_operating_point(
                    model, geometry, request, scales,
                    active_params, options
                );

            require_valid(
                active.solution, label + " incidence"
            );
            ++counts.incidence;

            const auto zero_diffusion =
                evaluate_diffusion_operating_point(
                    model, geometry, request, scales,
                    active_params, zero_diffusion_params,
                    options
                );

            require_valid(
                zero_diffusion.validation,
                label + " zero diffusion"
            );

            require_equal(
                active.solution,
                zero_diffusion.validation
            );

            ++counts.zero_diffusion;

            const auto combined =
                evaluate_diffusion_operating_point(
                    model, geometry, request, scales,
                    active_params, diffusion_params,
                    options
                );

            require_valid(
                combined.validation,
                label + " combined diffusion"
            );

            require(
                combined.independently_validated,
                label + ": independent diffusion validation failed"
            );

            require(
                close(
                    combined.total_entropy_generation,
                    combined.baseline_entropy_generation +
                    combined.incidence_entropy_generation +
                    combined.diffusion_entropy_generation
                ),
                label + ": entropy decomposition mismatch"
            );

            ++counts.combined;

            if (combined.diffusion.loss.specific_loss > 0.0)
                ++counts.active_diffusion;

            const auto baseline_map =
                evaluate_compressor_map_point(
                    model, geometry, baseline,
                    scales, reference
                );

            const auto active_map =
                evaluate_compressor_map_point(
                    model, geometry, active.solution,
                    scales, reference
                );

            const auto combined_map =
                evaluate_compressor_map_point(
                    model, geometry, combined.validation,
                    scales, reference
                );

            require(
                close(
                    baseline_map.corrected_mass_flow,
                    combined_map.corrected_mass_flow
                ) &&
                close(
                    baseline_map.corrected_speed_rpm,
                    combined_map.corrected_speed_rpm
                ),
                label + ": combined corrected coordinates changed"
            );

            ++counts.mapped;

            // Inlet total conditions are identical for
            // the baseline and incidence-aware solutions.
            require(
                close(
                    baseline_map.corrected_mass_flow,
                    active_map.corrected_mass_flow
                ) &&
                close(
                    baseline_map.corrected_speed_rpm,
                    active_map.corrected_speed_rpm
                ),
                label + ": corrected coordinates changed"
            );

            // At the calibrated reference point,
            // incidence loss must vanish.
            if (flow == 5.0 && speed == 12000.0) {
                require(
                    close(
                        active.incidence.loss.specific_loss,
                        0.0
                    ),
                    label + ": nonzero design incidence"
                );

                require_equal(
                    baseline, active.solution
                );
            }

            // The entropy model is additive.
            require(
                close(
                    active.total_entropy_generation,
                    active.baseline_entropy_generation +
                    active.incidence_entropy_generation
                ),
                label + ": entropy decomposition mismatch"
            );

            const auto& a = baseline_map;
            const auto& b = active_map;

            csv
                << name << ','
                << flow << ','
                << speed << ','
                << blade_angle << ','
                << active.incidence.loss.incidence_angle_deg << ','
                << inlet.value.rho << ','
                << active.incidence.inlet_relative_velocity << ','
                << active.incidence.loss.specific_loss << ','
                << active.baseline_entropy_generation << ','
                << active.incidence_entropy_generation << ','
                << active.total_entropy_generation << ','
                << baseline.outlet_temperature << ','
                << active.solution.outlet_temperature << ','
                << active.solution.outlet_temperature -
                   baseline.outlet_temperature << ','
                << baseline.outlet_pressure << ','
                << active.solution.outlet_pressure << ','
                << active.solution.outlet_pressure -
                   baseline.outlet_pressure << ','
                << a.static_pressure_ratio << ','
                << b.static_pressure_ratio << ','
                << b.static_pressure_ratio -
                   a.static_pressure_ratio << ','
                << a.total_pressure_ratio << ','
                << b.total_pressure_ratio << ','
                << b.total_pressure_ratio -
                   a.total_pressure_ratio << ','
                << a.specific_work << ','
                << b.specific_work << ','
                << b.specific_work -
                   a.specific_work << ','
                << a.shaft_power << ','
                << b.shaft_power << ','
                << b.shaft_power -
                   a.shaft_power << ','
                << a.torque << ','
                << b.torque << ','
                << b.torque -
                   a.torque << ','
                << a.corrected_mass_flow << ','
                << b.corrected_mass_flow << ','
                << a.corrected_speed_rpm << ','
                << b.corrected_speed_rpm << ','
                << a.full_residual_norm << ','
                << b.full_residual_norm << ','
                << combined.diffusion.loss.diffusion_ratio << ','
                << combined.diffusion.loss.specific_loss << ','
                << combined.diffusion_entropy_generation << ','
                << combined.total_entropy_generation << ','
                << combined_map.static_pressure_ratio << ','
                << combined_map.total_pressure_ratio << ','
                << combined.validation.outlet_temperature << ','
                << combined.validation.specific_work << ','
                << combined.validation.shaft_power << ','
                << combined.validation.validation_residual_norm
                << '\n';
        }
    }

    std::cout
        << name
        << ": points=" << counts.points
        << " baseline=" << counts.baseline
        << " zero=" << counts.zero
        << " incidence=" << counts.incidence
        << " mapped=" << counts.mapped
        << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-diffusion-map-comparison/"
                  "comparison.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path()
            );

        std::ofstream csv(output);

        require(
            static_cast<bool>(csv),
            "Cannot open benchmark CSV"
        );

        csv << std::setprecision(17);
        write_header(csv);

        const thermomodel::IdealGas air_like;

        thermomodel::PRParameters pr_params;

        thermomodel::IdealGasParameters ideal_params;
        ideal_params.R = Ru / pr_params.molar_mass;
        ideal_params.cp = methane_cp;
        ideal_params.Tref = pr_params.Tref;
        ideal_params.pref = pr_params.pref;

        const thermomodel::IdealGas methane_ideal(
            ideal_params
        );

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                methane_cp
            );

        const thermomodel::PengRobinson methane_pr(
            caloric, pr_params
        );

        Counts air_counts;
        Counts ideal_counts;
        Counts pr_counts;

        sweep(
            air_like, "AirLikeIdeal",
            csv, air_counts
        );

        sweep(
            methane_ideal, "MethaneIdeal",
            csv, ideal_counts
        );

        sweep(
            methane_pr, "MethanePR",
            csv, pr_counts
        );

        csv.close();

        const auto accepted = [](const Counts& c) {
            return c.points == 35 &&
                   c.baseline == 35 &&
                   c.zero == 35 &&
                   c.incidence == 35 &&
                   c.mapped == 35 &&
                   c.combined == 35 &&
                   c.zero_diffusion == 35;
        };

        require(
            accepted(air_counts) &&
            accepted(ideal_counts) &&
            accepted(pr_counts),
            "Incidence-map benchmark acceptance failed"
        );

        std::cout
            << "PASS: 105 paired operating points\n"
            << "PASS: 315 primary nonlinear solutions\n"
            << "PASS: 105 zero-diffusion regressions\n"
            << "PASS: 105 zero-coefficient regressions\n"
            << "CSV: " << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
