#include <turbocompgpu/complete_stage_friction.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    const char* label,
    double actual,
    double expected,
    double tolerance = 1.0e-8)
{
    const double scale = std::max({
        1.0, std::abs(actual), std::abs(expected)
    });

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        std::cerr << label
                  << " actual=" << actual
                  << " expected=" << expected << '\n';
        throw std::runtime_error(label);
    }
}

template<class Model>
void run_case(const Model& model, const char* name)
{
    const ImpellerSystemGeometry geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearResidualScales scales{
        50000.0, 100.0, 500000.0
    };

    const OperatingPointRequest request{
        {
            300.0,
            500000.0,
            5.0,
            12000.0,
            0.90,
            20.0,
            20.0
        },
        {
            340.0,
            std::log(700000.0 / scales.pressure)
        }
    };

    NonlinearSolverOptions solver_options;
    solver_options.residual_tolerance = 1.0e-13;

    // Match the impeller's effective exit area exactly.
    const double r2 = geometry.exit.diameter / 2.0;

    const double area2 = effective_flow_area(
        impeller_exit_area(geometry.exit),
        geometry.exit.blockage);

    const double width = area2 /
        (2.0 * std::numbers::pi * r2);

    const VanelessDiffuserGeometry diffuser_geometry{
        r2,
        1.5 * r2,
        width,
        0.0
    };

    FrictionalDiffuserOptions diffuser_options;
    diffuser_options.darcy_friction_factor = 0.01;
    diffuser_options.radial_steps = 512;

    const auto result = evaluate_frictional_complete_stage(
        model,
        model,
        geometry,
        diffuser_geometry,
        request,
        scales,
        solver_options,
        diffuser_options);

    if (!result.complete())
        std::cerr << name << ": " << result.message << '\n';

    require(result.complete(), "Complete-stage solve failed");
    require(result.impeller_valid(), "Impeller validation failed");
    require(result.diffuser.has_value(), "Missing diffuser");
    require(result.station1_stagnation.has_value(),
            "Missing station 1 stagnation");
    require(result.station3_stagnation.has_value(),
            "Missing station 3 stagnation");

    const auto& d = *result.diffuser;

    close("mass conservation", d.mass_residual, 0.0);
    close("diffuser energy conservation",
          d.stagnation_enthalpy_residual, 0.0);
    close("combined energy conservation",
          result.combined_energy_residual, 0.0,
          1.0e-7);
    close("entropy balance",
          d.entropy_balance_residual, 0.0,
          1.0e-7);
    close("angular momentum balance",
          d.angular_momentum_balance_residual, 0.0,
          1.0e-7);

    require(d.entropy_change >= -1.0e-8,
            "Negative diffuser entropy production");

    close("fluid power",
          result.fluid_power,
          request.operating.mass_flow *
              result.fluid_specific_work);

    // The diffuser must not add shaft work.
    close("stage work",
          result.fluid_specific_work,
          result.impeller.specific_work,
          1.0e-7);

    // Explicit geometry mismatch must fail without
    // discarding the validated impeller result.
    auto invalid_geometry = diffuser_geometry;
    invalid_geometry.outlet_width *= 1.05;

    const auto invalid = evaluate_frictional_complete_stage(
        model,
        model,
        geometry,
        invalid_geometry,
        request,
        scales,
        solver_options,
        diffuser_options);

    require(
        invalid.status == CompleteStageStatus::DiffuserFailed,
        "Geometry mismatch was not rejected");

    require(invalid.impeller_valid(),
            "Failed diffuser invalidated impeller");

    std::cout
        << name
        << " p3=" << d.outlet_static.p
        << " T3=" << d.outlet_static.T
        << " total_pressure_ratio="
        << result.total_to_total_pressure_ratio
        << " fluid_power=" << result.fluid_power
        << " entropy_change=" << d.entropy_change
        << " energy_residual="
        << result.combined_energy_residual
        << '\n';
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas ideal;
        run_case(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(caloric, parameters);
        run_case(pr, "PengRobinson");

        std::cout << "Frictional complete-stage integration: PASS\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
