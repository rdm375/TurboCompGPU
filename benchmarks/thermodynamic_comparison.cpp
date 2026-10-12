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

struct Counts {
    int attempted = 0;
    int converged = 0;
    int validated = 0;
    int physical = 0;
    int mapped = 0;
};

void write_header(std::ofstream& out)
{
    out
        << "model,mass_flow,rotational_speed_rpm,"
        << "converged,validated,basic_physics_valid,"
        << "iterations,residual_norm,"
        << "inlet_temperature,inlet_pressure,"
        << "inlet_density,inlet_enthalpy,inlet_entropy,"
        << "inlet_cp,inlet_Z,"
        << "outlet_temperature,outlet_pressure,"
        << "outlet_density,outlet_enthalpy,outlet_entropy,"
        << "outlet_cp,outlet_Z,"
        << "inlet_total_temperature,inlet_total_pressure,"
        << "outlet_total_temperature,outlet_total_pressure,"
        << "corrected_mass_flow,corrected_speed_rpm,"
        << "static_pressure_ratio,total_pressure_ratio,"
        << "specific_work,shaft_power,torque,"
        << "delta_static_enthalpy,delta_static_entropy,"
        << "entropy_generation,"
        << "full_residual_norm,"
        << "inlet_stagnation_h_residual,"
        << "inlet_stagnation_s_residual,"
        << "outlet_stagnation_h_residual,"
        << "outlet_stagnation_s_residual\n";
}

void sweep(
    const thermomodel::Model& model,
    const std::string& name,
    std::ofstream& out,
    Counts& counts)
{
    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    constexpr std::array<double, 7> flows {
        2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
    };

    constexpr std::array<double, 5> speeds {
        9000.0, 10500.0, 12000.0, 13500.0, 15000.0
    };

    for (double speed : speeds) {
        for (double flow : flows) {
            const NonlinearOperatingPoint operating {
                300.0, 500000.0,
                flow, speed,
                0.90, 20.0, 20.0
            };

            const NonlinearUnknowns guess {
                340.0,
                std::log(700000.0 / scales.pressure)
            };

            ++counts.attempted;

            const auto point = evaluate_operating_point(
                model,
                geometry,
                {operating, guess},
                scales,
                options
            );

            counts.converged += point.converged;
            counts.validated += point.independently_validated;
            counts.physical += point.basic_physics_valid;

            if (!point.converged ||
                !point.independently_validated ||
                !point.basic_physics_valid) {
                std::cerr
                    << "FAILED: " << name
                    << " flow=" << flow
                    << " speed=" << speed
                    << " message=" << point.message
                    << '\n';
                continue;
            }

            const auto map = evaluate_compressor_map_point(
                model, geometry, point, scales, reference
            );

            ++counts.mapped;

            const auto inlet = model.evaluate_Tp(
                operating.inlet_temperature,
                operating.inlet_pressure
            );

            const auto outlet = model.evaluate_Tp(
                point.outlet_temperature,
                point.outlet_pressure
            );

            if (!inlet || !outlet)
                throw std::runtime_error(
                    "Thermodynamic state reconstruction failed"
                );

            const auto& a = inlet.value;
            const auto& b = outlet.value;

            out
                << name << ','
                << flow << ','
                << speed << ','
                << point.converged << ','
                << point.independently_validated << ','
                << point.basic_physics_valid << ','
                << point.iterations << ','
                << point.residual_norm << ','
                << a.T << ',' << a.p << ','
                << a.rho << ',' << a.h << ',' << a.s << ','
                << a.cp << ',' << a.Z << ','
                << b.T << ',' << b.p << ','
                << b.rho << ',' << b.h << ',' << b.s << ','
                << b.cp << ',' << b.Z << ','
                << map.inlet_total_temperature << ','
                << map.inlet_total_pressure << ','
                << map.outlet_total_temperature << ','
                << map.outlet_total_pressure << ','
                << map.corrected_mass_flow << ','
                << map.corrected_speed_rpm << ','
                << map.static_pressure_ratio << ','
                << map.total_pressure_ratio << ','
                << map.specific_work << ','
                << map.shaft_power << ','
                << map.torque << ','
                << b.h - a.h << ','
                << b.s - a.s << ','
                << map.entropy_generation << ','
                << map.full_residual_norm << ','
                << map.inlet_stagnation_enthalpy_residual << ','
                << map.inlet_stagnation_entropy_residual << ','
                << map.outlet_stagnation_enthalpy_residual << ','
                << map.outlet_stagnation_entropy_residual
                << '\n';
        }
    }

    std::cout
        << name
        << ": attempted=" << counts.attempted
        << " converged=" << counts.converged
        << " validated=" << counts.validated
        << " physical=" << counts.physical
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
                : "benchmarks/results-thermodynamic-comparison/"
                  "comparison.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path()
            );

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open CSV output");

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

        std::cout << std::setprecision(12)
                  << "Matched methane R = "
                  << ideal_params.R << " J/(kg K)\n"
                  << "Matched ideal-gas cp = "
                  << methane_cp << " J/(kg K)\n";

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
            return c.attempted == 35 &&
                   c.converged == 35 &&
                   c.validated == 35 &&
                   c.physical == 35 &&
                   c.mapped == 35;
        };

        if (!accepted(air_counts) ||
            !accepted(ideal_counts) ||
            !accepted(pr_counts)) {
            throw std::runtime_error(
                "Thermodynamic comparison acceptance failed"
            );
        }

        std::cout
            << "PASS: 105/105 operating points\n"
            << "CSV: " << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
