#include <turbocompgpu/operating_diffusion_boundary.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <algorithm>
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

constexpr std::array<double, 5> speeds {
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

const ImpellerSystemGeometry geometry {
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const NonlinearResidualScales scales {
    50000.0, 100.0, 500000.0
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

bool close(double a, double b, double tolerance = 1e-8)
{
    return std::isfinite(a) &&
           std::isfinite(b) &&
           std::abs(a - b) <= tolerance *
               std::max({1.0, std::abs(a), std::abs(b)});
}

OperatingPointRequest make_request(double flow, double speed)
{
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

    return {operating, guess};
}

struct Counts {
    int cases = 0;
    int converged = 0;
    int unbracketed = 0;
    int failed = 0;
};

void sweep(
    const thermomodel::Model& model,
    const std::string& name,
    std::ofstream& csv,
    Counts& counts)
{
    NonlinearSolverOptions solver_options;
    solver_options.residual_tolerance = 1.0e-13;

    DiffusionBoundaryOptions boundary_options;
    boundary_options.flow_tolerance = 1.0e-6;
    boundary_options.diffusion_tolerance = 1.0e-8;
    boundary_options.max_iterations = 80;

    const auto inlet = model.evaluate_Tp(300.0, 500000.0);

    require(
        static_cast<bool>(inlet),
        name + ": inlet evaluation failed"
    );

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

    const IncidenceLossParameters incidence_params {
        blade_angle, 1.0
    };

    const DiffusionLossParameters diffusion_params {
        1.0, 1.5
    };

    std::cout
        << "\n" << name
        << " blade_angle_deg=" << blade_angle
        << '\n';

    for (double speed : speeds) {
        ++counts.cases;

        const auto evaluate = [&](double flow) {
            const auto point =
                evaluate_operating_diffusion_boundary(
                    model,
                    geometry,
                    make_request(flow, speed),
                    scales,
                    incidence_params,
                    diffusion_params,
                    solver_options
                );

            return point.diffusion_ratio;
        };

        const auto boundary = solve_diffusion_boundary(
            evaluate,
            2.0,
            8.0,
            diffusion_params.critical_diffusion_ratio,
            boundary_options
        );

        if (boundary.status ==
            DiffusionBoundaryStatus::not_bracketed) {

            ++counts.unbracketed;

            const double ratio_low = evaluate(2.0);
            const double ratio_high = evaluate(8.0);

            csv
                << name << ','
                << speed << ','
                << blade_angle << ','
                << "not_bracketed,"
                << ",,,,"
                << boundary.evaluations << ','
                << boundary.iterations << ','
                << ratio_low << ','
                << ratio_high
                << ",,,,,,,\n";

            std::cout
                << "  N=" << speed
                << " unbracketed"
                << " Dr(2)=" << ratio_low
                << " Dr(8)=" << ratio_high
                << '\n';

            continue;
        }

        if (!boundary.converged()) {
            ++counts.failed;

            throw std::runtime_error(
                name + ": boundary solver did not converge"
            );
        }

        const double flow = boundary.critical_mass_flow;

        const auto point =
            evaluate_operating_diffusion_boundary(
                model,
                geometry,
                make_request(flow, speed),
                scales,
                incidence_params,
                diffusion_params,
                solver_options
            );

        const auto combined =
            evaluate_diffusion_operating_point(
                model,
                geometry,
                make_request(flow, speed),
                scales,
                incidence_params,
                diffusion_params,
                solver_options
            );

        require(
            combined.independently_validated &&
            combined.validation.converged &&
            combined.validation.basic_physics_valid,
            name + ": coupled boundary validation failed"
        );

        require(
            std::abs(point.diffusion_ratio - 1.5) <=
                boundary_options.diffusion_tolerance,
            name + ": diffusion ratio mismatch"
        );

        require(
            boundary.upper_mass_flow -
                boundary.lower_mass_flow <=
                boundary_options.flow_tolerance,
            name + ": bracket tolerance failed"
        );

        // At the threshold the loss is mathematically zero.
        // The independently solved states should agree.
        require(
            close(
                point.incidence_point.solution.outlet_pressure,
                combined.validation.outlet_pressure
            ),
            name + ": pressure mismatch at activation boundary"
        );

        require(
            close(
                point.incidence_point.solution.outlet_temperature,
                combined.validation.outlet_temperature
            ),
            name + ": temperature mismatch at activation boundary"
        );

        // Verify that the root separates active and
        // inactive diffusion regimes.
        //
        // Use a perturbation much larger than the root
        // tolerance, but small relative to the operating
        // mass flow.
        const double perturbation = 1.0e-3;

        const double ratio_below =
            evaluate(flow - perturbation);

        const double ratio_above =
            evaluate(flow + perturbation);

        require(
            ratio_below >
                diffusion_params.critical_diffusion_ratio,
            name + ": diffusion not active below boundary"
        );

        require(
            ratio_above <
                diffusion_params.critical_diffusion_ratio,
            name + ": diffusion not inactive above boundary"
        );

        const double U1 =
            point.incidence_point.incidence.inlet_blade_speed;

        const double phi1 =
            point.inlet_meridional_velocity / U1;

        const double pressure_difference =
            combined.validation.outlet_pressure -
            point.incidence_point.solution.outlet_pressure;

        const double temperature_difference =
            combined.validation.outlet_temperature -
            point.incidence_point.solution.outlet_temperature;

        ++counts.converged;

        csv
            << name << ','
            << speed << ','
            << blade_angle << ','
            << "converged,"
            << flow << ','
            << point.diffusion_ratio << ','
            << boundary.diffusion_residual << ','
            << boundary.upper_mass_flow -
               boundary.lower_mass_flow << ','
            << boundary.evaluations << ','
            << boundary.iterations << ','
            << ",,"
            << phi1 << ','
            << inlet.value.rho << ','
            << point.specific_loss << ','
            << combined.diffusion.loss.specific_loss << ','
            << pressure_difference << ','
            << temperature_difference << ','
            << combined.validation.validation_residual_norm
            << '\n';

        std::cout
            << "  N=" << speed
            << " mcrit=" << flow
            << " Dr=" << point.diffusion_ratio
            << " phi1=" << phi1
            << " evaluations=" << boundary.evaluations
            << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-diffusion-boundary/"
                  "boundary.csv";

        if (!output.parent_path().empty())
            std::filesystem::create_directories(
                output.parent_path()
            );

        std::ofstream csv(output);

        require(
            static_cast<bool>(csv),
            "Cannot open boundary CSV"
        );

        csv << std::setprecision(17);

        csv
            << "model,rotational_speed_rpm,blade_angle_deg,"
            << "status,critical_mass_flow,diffusion_ratio,"
            << "diffusion_residual,bracket_width,"
            << "evaluations,iterations,"
            << "diffusion_ratio_at_2,diffusion_ratio_at_8,"
            << "phi1,inlet_density,"
            << "incidence_diffusion_loss,"
            << "coupled_diffusion_loss,"
            << "delta_outlet_pressure,"
            << "delta_outlet_temperature,"
            << "coupled_full_residual_norm\n";

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

        Counts counts;

        sweep(air_like, "AirLikeIdeal", csv, counts);
        sweep(methane_ideal, "MethaneIdeal", csv, counts);
        sweep(methane_pr, "MethanePR", csv, counts);

        csv.flush();

        require(
            static_cast<bool>(csv),
            "Failed to write boundary CSV"
        );

        std::cout
            << "\nBoundary summary:\n"
            << "  cases=" << counts.cases << '\n'
            << "  converged=" << counts.converged << '\n'
            << "  unbracketed=" << counts.unbracketed << '\n'
            << "  failed=" << counts.failed << '\n';

        require(
            counts.cases == 15 &&
            counts.failed == 0 &&
            counts.converged + counts.unbracketed == 15,
            "Incomplete boundary sweep"
        );

        std::cout << "PASS: diffusion boundary sweep\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
