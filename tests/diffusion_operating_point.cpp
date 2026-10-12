#include <turbocompgpu/diffusion_operating_point.hpp>

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

void close(
    double actual,
    double expected,
    double tolerance = 1e-9)
{
    const double scale = std::max({
        1.0,
        std::abs(actual),
        std::abs(expected)
    });

    check(
        std::isfinite(actual) &&
        std::isfinite(expected) &&
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
        50000.0, 100.0, 500000.0
    };

    const auto inlet = model.evaluate_Tp(
        300.0, 500000.0);

    check(
        static_cast<bool>(inlet),
        "Inlet evaluation failed");

    const OperatingIncidenceConditions design {
        5.0, 12000.0, inlet.value.rho, 20.0
    };

    const double blade_angle =
        calibrate_inlet_blade_angle_deg(
            geometry, design);

    const NonlinearUnknowns guess {
        340.0,
        std::log(700000.0 / scales.pressure)
    };

    const OperatingPointRequest request {
        {
            300.0,
            500000.0,
            5.0,
            12000.0,
            0.90,
            20.0,
            20.0
        },
        guess
    };

    // Reference: previously validated incidence solver.
    const auto reference =
        evaluate_incidence_operating_point(
            model,
            geometry,
            request,
            scales,
            {blade_angle, 1.0});

    check(
        reference.solution.converged &&
        reference.solution.independently_validated,
        "Incidence reference failed");

    // Zero diffusion must recover incidence-only behavior.
    const auto zero =
        evaluate_diffusion_operating_point(
            model,
            geometry,
            request,
            scales,
            {blade_angle, 1.0},
            {0.0, 1.5});

    check(
        zero.independently_validated &&
        zero.validation.basic_physics_valid,
        "Zero-diffusion validation failed");

    close(
        zero.diffusion_entropy_generation,
        0.0);

    close(
        zero.validation.outlet_temperature,
        reference.solution.outlet_temperature);

    close(
        zero.validation.outlet_pressure,
        reference.solution.outlet_pressure);

    close(
        zero.validation.specific_work,
        reference.solution.specific_work);

    // Activate a low threshold so the new coupling
    // is exercised at the design operating point.
    const auto active =
        evaluate_diffusion_operating_point(
            model,
            geometry,
            request,
            scales,
            {blade_angle, 1.0},
            {0.01, 1.0});

    check(
        active.independently_validated &&
        active.validation.basic_physics_valid,
        "Active diffusion validation failed");

    close(
        active.total_entropy_generation,
        active.baseline_entropy_generation +
        active.incidence_entropy_generation +
        active.diffusion_entropy_generation);

    close(
        active.diffusion_entropy_generation,
        active.diffusion.loss.specific_loss /
        active.inlet_stagnation_temperature);

    check(
        active.coupled_residual_norm <= 1e-8,
        "Coupled residual too large");

    // Reconstruct relative outlet velocity independently
    // from the converged EOS density and slip relation.
    const double pi = std::acos(-1.0);

    const double area =
        pi * geometry.exit.diameter *
        geometry.exit.width *
        (1.0 - geometry.exit.blockage);

    const double vm2 =
        request.operating.mass_flow /
        (active.validation.outlet_density * area);

    const double omega =
        2.0 * pi *
        request.operating.rotational_speed_rpm / 60.0;

    const double U2 =
        omega * geometry.impeller.exit_diameter / 2.0;

    const double beta =
        geometry.impeller.exit_blade_angle_deg *
        pi / 180.0;

    const double whirl =
        request.operating.slip_factor *
        (U2 - vm2 / std::tan(beta));

    const double W2 =
        std::hypot(vm2, whirl - U2);

    close(
        active.diffusion.loss.outlet_relative_velocity,
        W2);

    std::cout
        << name
        << " ratio=" << active.diffusion.loss.diffusion_ratio
        << " diffusion_loss="
        << active.diffusion.loss.specific_loss
        << " diffusion_entropy="
        << active.diffusion_entropy_generation
        << " residual="
        << active.coupled_residual_norm
        << '\n';
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
            ideal_params);

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane_pr(
            caloric, pr_params);

        check_model(air_like, "AirLikeIdeal");
        check_model(methane_ideal, "MethaneIdeal");
        check_model(methane_pr, "MethanePR");

        std::cout
            << "PASS: coupled diffusion operating points\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
