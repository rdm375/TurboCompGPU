#include <turbocompgpu/incidence_operating_point.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

void check(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void check_close(
    double actual,
    double expected,
    double tolerance = 1e-10)
{
    const double scale = std::max(
        {1.0, std::abs(actual), std::abs(expected)}
    );

    check(
        std::isfinite(actual) &&
        std::abs(actual - expected) <= tolerance * scale,
        "Numerical mismatch"
    );
}

void check_model(
    const thermomodel::Model& model,
    const char* name)
{
    const ImpellerSystemGeometry geometry {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearResidualScales scales {
        50000.0,
        100.0,
        500000.0
    };

    const auto inlet = model.evaluate_Tp(
        300.0, 500000.0
    );

    check(
        static_cast<bool>(inlet),
        "Inlet thermodynamic evaluation failed"
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
        blade_angle,
        1.0
    };

    const NonlinearUnknowns guess {
        340.0,
        std::log(700000.0 / scales.pressure)
    };

    const auto make_request = [&](double mass_flow) {
        return OperatingPointRequest {
            {
                300.0,
                500000.0,
                mass_flow,
                12000.0,
                0.90,
                20.0,
                20.0
            },
            guess
        };
    };

    for (double mass_flow : {3.0, 5.0, 7.0}) {
        const auto request = make_request(mass_flow);

        // Reference solution from the existing solver.
        const auto baseline = evaluate_operating_point(
            model, geometry, request, scales
        );

        check(
            baseline.converged &&
            baseline.independently_validated &&
            baseline.basic_physics_valid,
            "Baseline solve failed"
        );

        // Zero incidence coefficient must reproduce
        // the original nonlinear problem.
        const auto zero = evaluate_incidence_operating_point(
            model,
            geometry,
            request,
            scales,
            {blade_angle, 0.0}
        );

        check(
            zero.solution.converged &&
            zero.solution.independently_validated &&
            zero.solution.basic_physics_valid,
            "Zero-loss solve failed"
        );

        check_close(
            zero.incidence_entropy_generation,
            0.0
        );

        check_close(
            zero.solution.outlet_temperature,
            baseline.outlet_temperature
        );

        check_close(
            zero.solution.outlet_pressure,
            baseline.outlet_pressure
        );

        check_close(
            zero.solution.specific_work,
            baseline.specific_work
        );

        // Activate incidence loss.
        const auto with_incidence =
            evaluate_incidence_operating_point(
                model,
                geometry,
                request,
                scales,
                incidence_params
            );

        check(
            with_incidence.solution.converged &&
            with_incidence.solution.independently_validated &&
            with_incidence.solution.basic_physics_valid,
            "Incidence-aware solve failed"
        );

        check_close(
            with_incidence.total_entropy_generation,
            with_incidence.baseline_entropy_generation +
                with_incidence.incidence_entropy_generation
        );

        check_close(
            with_incidence.solution.entropy_generation,
            with_incidence.total_entropy_generation
        );

        if (mass_flow == 5.0) {
            // Design condition: no incidence penalty.
            check_close(
                with_incidence.incidence_entropy_generation,
                0.0
            );

            check_close(
                with_incidence.solution.outlet_pressure,
                baseline.outlet_pressure
            );
        }
        else {
            check(
                with_incidence.incidence_entropy_generation > 0.0,
                "Off-design entropy contribution not positive"
            );
        }

        std::cout
            << name
            << " flow=" << mass_flow
            << " incidence_deg="
            << with_incidence.incidence.loss.incidence_angle_deg
            << " incidence_loss_J_kg="
            << with_incidence.incidence.loss.specific_loss
            << " incidence_entropy_J_kgK="
            << with_incidence.incidence_entropy_generation
            << " pressure_ratio="
            << with_incidence.solution.pressure_ratio
            << '\n';
    }
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas air_like;

        thermomodel::PRParameters pr_params;

        thermomodel::IdealGasParameters ideal_params;
        ideal_params.R =
            8.31446261815324 / pr_params.molar_mass;
        ideal_params.cp = 2200.0;
        ideal_params.Tref = pr_params.Tref;
        ideal_params.pref = pr_params.pref;

        const thermomodel::IdealGas methane_ideal(
            ideal_params
        );

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0
            );

        const thermomodel::PengRobinson methane_pr(
            caloric, pr_params
        );

        check_model(air_like, "AirLikeIdeal");
        check_model(methane_ideal, "MethaneIdeal");
        check_model(methane_pr, "MethanePR");

        std::cout
            << "PASS: incidence-aware operating points\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
