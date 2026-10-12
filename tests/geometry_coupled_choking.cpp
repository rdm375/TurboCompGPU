#include <turbocompgpu/geometry_coupled_choking.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double tolerance,
    const char* message)
{
    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <=
            tolerance *
            std::max(1.0, std::abs(expected)),
        message);
}

using namespace turbocompgpu;

const ImpellerSystemGeometry impeller_geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

// Deliberately select a narrow local passage height
// to give an effective area near the Stage 5 reference
// value of 0.002 m^2.
//
// These are illustrative dimensions, not a calibrated
// compressor geometry.
const ThroatGeometry throat_geometry{
    12,
    0.15,
    0.0035,
    45.0,
    0.002,
    0.10
};

template<class Model>
void run_case(
    const Model& model,
    const char* name,
    double rpm)
{
    const GeometryCoupledChokingInputs input{
        300.0,
        500000.0,
        rpm,
        20.0,
        throat_geometry
    };

    const auto result =
        evaluate_geometry_coupled_choking(
            model, model, impeller_geometry, input);

    const auto& coupled = result.coupled;
    const auto& critical = coupled.throat.critical;

    require(
        coupled.critical_mass_flow > 0.0,
        "Nonpositive critical mass flow");

    require(
        critical.density > 0.0 &&
        critical.sound_speed > 0.0,
        "Invalid critical thermodynamic state");

    // Independent geometric expression.
    const double beta =
        throat_geometry.blade_angle_deg *
        std::numbers::pi / 180.0;

    const double pitch =
        2.0 * std::numbers::pi *
        throat_geometry.radius /
        throat_geometry.blade_count;

    const double expected_area =
        throat_geometry.blade_count *
        throat_geometry.passage_height *
        (pitch * std::sin(beta) -
         throat_geometry.blade_thickness) *
        (1.0 - throat_geometry.blockage);

    close(
        result.geometry.effective_area,
        expected_area,
        1e-12,
        "Geometry-derived area mismatch");

    // Verify sonic flow in throat-local coordinates.
    close(
        result.throat_local_velocity.normal,
        critical.sound_speed,
        1e-12,
        "Throat-normal sonic condition");

    close(
        result.throat_local_velocity.tangential,
        0.0,
        1e-12,
        "Blade-aligned throat tangential velocity");

    // Independently verify the impeller-coordinate
    // velocity components.
    close(
        result.throat_relative_velocity.meridional,
        critical.sound_speed * std::sin(beta),
        1e-12,
        "Meridional sonic component");

    close(
        result.throat_relative_velocity.circumferential,
        -critical.sound_speed * std::cos(beta),
        1e-12,
        "Circumferential sonic component");

    // Full relative kinetic energy, without adding
    // the circumferential component a second time.
    const double W2 =
        result.throat_relative_velocity.meridional *
            result.throat_relative_velocity.meridional +
        result.throat_relative_velocity.circumferential *
            result.throat_relative_velocity.circumferential;

    close(
        W2,
        critical.sound_speed * critical.sound_speed,
        1e-10,
        "Relative kinetic-energy closure");

    // Stage 6C energy-residual decomposition.
    const auto effective = model.evaluate_Tp(
        coupled.effective_stagnation_temperature,
        coupled.effective_stagnation_pressure);

    require(
        static_cast<bool>(effective),
        "Effective stagnation EOS evaluation failed");

    const double sonic_energy_residual =
        critical.enthalpy +
        0.5 * critical.sound_speed *
            critical.sound_speed -
        effective.value.h;

    const double effective_state_residual =
        effective.value.h -
        0.5 * result.throat_blade_speed *
            result.throat_blade_speed -
        coupled.inlet_rothalpy;

    const double reconstructed_residual =
        sonic_energy_residual +
        effective_state_residual;

    std::cout
        << std::setprecision(15)
        << name
        << " rpm=" << rpm
        << " energy diagnostics:\n"
        << "  sonic energy residual     = "
        << sonic_energy_residual << " J/kg\n"
        << "  effective state residual  = "
        << effective_state_residual << " J/kg\n"
        << "  reconstructed residual    = "
        << reconstructed_residual << " J/kg\n"
        << "  direct rothalpy residual  = "
        << result.rothalpy_residual << " J/kg\n"
        << "  sonic solver residual     = "
        << critical.energy_residual << " J/kg\n";

    const double energy_scale =
        std::max(
            1.0,
            std::abs(coupled.inlet_rothalpy));

    require(
        std::abs(result.rothalpy_residual) <
            1e-8 * energy_scale,
        "Rothalpy conservation residual too large");

    require(
        std::abs(result.sonic_residual) <
            1e-10,
        "Sonic residual too large");

    require(
        std::abs(result.mass_flow_residual) <
            1e-8 *
            std::max(1.0, coupled.critical_mass_flow),
        "Throat mass-flow residual too large");

    require(
        std::abs(result.continuity_residual) <
            1e-8 *
            std::max(1.0, coupled.critical_mass_flow),
        "Inlet continuity residual too large");

    require(
        std::abs(coupled.mass_flow_residual) <
            1e-8 *
            std::max(1.0, coupled.critical_mass_flow),
        "Coupled root residual too large");

    // Independent thermodynamic state checks.
    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    const auto throat = model.evaluate_Tp(
        critical.temperature,
        critical.pressure);

    require(
        static_cast<bool>(inlet) &&
        static_cast<bool>(throat),
        "EOS state evaluation failed");

    require(
        std::abs(throat.value.s - inlet.value.s) <
            1e-8,
        "Isentropic throat entropy mismatch");

    const double direct_rothalpy =
        throat.value.h +
        0.5 * W2 -
        0.5 * result.throat_blade_speed *
            result.throat_blade_speed;

    require(
        std::abs(
            direct_rothalpy -
            coupled.inlet_rothalpy) <
            1e-8 * energy_scale,
        "Independent EOS rothalpy mismatch");

    std::cout
        << std::setprecision(12)
        << name
        << " rpm=" << rpm
        << " mdot*=" << coupled.critical_mass_flow
        << " kg/s"
        << " area=" << result.geometry.effective_area
        << " m^2"
        << " Wm*="
        << result.throat_relative_velocity.meridional
        << " m/s"
        << " Wtheta*="
        << result.throat_relative_velocity.circumferential
        << " m/s"
        << " dI=" << result.rothalpy_residual
        << " J/kg"
        << " dmdot=" << result.mass_flow_residual
        << " kg/s\n";
}

void test_geometry_scaling()
{
    const auto base =
        evaluate_throat_geometry(throat_geometry);

    auto larger = throat_geometry;
    larger.passage_height *= 2.0;

    const auto doubled =
        evaluate_throat_geometry(larger);

    close(
        doubled.effective_area,
        2.0 * base.effective_area,
        1e-12,
        "Passage-height area scaling");

    std::cout
        << "PASS: throat area scales with passage height\n";
}

void test_stage5_equivalence()
{
    thermomodel::IdealGas model;

    const double rpm = 9000.0;

    const GeometryCoupledChokingInputs input{
        300.0,
        500000.0,
        rpm,
        20.0,
        throat_geometry
    };

    const auto integrated =
        evaluate_geometry_coupled_choking(
            model, model, impeller_geometry, input);

    const InletCoupledChokingInputs reference_input{
        input.inlet_temperature,
        input.inlet_pressure,
        input.rotational_speed_rpm,
        input.inlet_prewhirl,
        throat_geometry.radius,
        0.0,
        integrated.geometry.effective_area
    };

    const auto reference =
        evaluate_inlet_coupled_choking(
            model,
            model,
            impeller_geometry,
            reference_input);

    close(
        integrated.coupled.critical_mass_flow,
        reference.critical_mass_flow,
        1e-12,
        "Stage 5 zero-tangential-velocity equivalence");

    close(
        integrated.coupled.inlet_rothalpy,
        reference.inlet_rothalpy,
        1e-12,
        "Stage 5 rothalpy equivalence");

    std::cout
        << "PASS: Stage 5 zero-local-tangential equivalence\n";
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas air;

        thermomodel::PRParameters pr_params;
        pr_params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::IdealGasParameters methane_params;
        methane_params.R =
            8.31446261815324 / pr_params.molar_mass;
        methane_params.cp = 2200.0;
        methane_params.Tref = pr_params.Tref;
        methane_params.pref = pr_params.pref;

        thermomodel::IdealGas methane_ideal(
            methane_params);

        auto caloric =
            std::make_shared<
                thermomodel::ConstantCp>(2200.0);

        thermomodel::PengRobinson methane_pr(
            caloric, pr_params);

        test_geometry_scaling();
        test_stage5_equivalence();

        for (double rpm : {
                 0.0, 9000.0, 12000.0, 15000.0
             }) {
            run_case(air, "AirLikeIdeal", rpm);
            run_case(
                methane_ideal, "MethaneIdeal", rpm);
            run_case(
                methane_pr, "MethanePR", rpm);
        }

        std::cout
            << "PASS: geometry-coupled choking tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
