#include <turbocompgpu/incidence_coupled_choking.hpp>
#include <turbocompgpu/geometry_coupled_choking.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <limits>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double relative_tolerance,
    const char* message)
{
    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <=
            relative_tolerance *
            std::max(1.0, std::abs(expected)),
        message);
}

const ImpellerSystemGeometry impeller_geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

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
    double rpm,
    double Wtau)
{
    const IncidenceCoupledChokingInputs input{
        300.0,
        500000.0,
        rpm,
        20.0,
        throat_geometry,
        Wtau
    };

    const auto result =
        evaluate_incidence_coupled_choking(
            model, model, impeller_geometry, input);

    const auto& coupled = result.coupled;
    const auto& critical = coupled.throat.critical;

    const double beta =
        throat_geometry.blade_angle_deg *
        std::numbers::pi / 180.0;

    const double a = critical.sound_speed;

    require(
        coupled.critical_mass_flow > 0.0,
        "Nonpositive critical mass flow");

    close(
        result.throat_local_velocity.normal,
        a,
        1e-12,
        "Normal sonic condition");

    close(
        result.throat_local_velocity.tangential,
        Wtau,
        1e-12,
        "Prescribed tangential velocity");

    // Independent coordinate expressions.
    close(
        result.throat_relative_velocity.meridional,
        a * std::sin(beta) +
            Wtau * std::cos(beta),
        1e-12,
        "Meridional velocity transformation");

    close(
        result.throat_relative_velocity.circumferential,
        -a * std::cos(beta) +
            Wtau * std::sin(beta),
        1e-12,
        "Circumferential velocity transformation");

    const double W2 =
        result.throat_relative_velocity.meridional *
            result.throat_relative_velocity.meridional +
        result.throat_relative_velocity.circumferential *
            result.throat_relative_velocity.circumferential;

    close(
        W2,
        a * a + Wtau * Wtau,
        1e-12,
        "Kinetic-energy invariance");

    const double energy_scale =
        std::max(
            1.0,
            std::abs(coupled.inlet_rothalpy));

    require(
        std::abs(result.rothalpy_residual) <
            1e-8 * energy_scale,
        "Rothalpy conservation residual");

    require(
        std::abs(result.sonic_residual) < 1e-10,
        "Sonic residual");

    const double flow_scale =
        std::max(
            1.0,
            coupled.critical_mass_flow);

    require(
        std::abs(result.mass_flow_residual) <
            1e-8 * flow_scale,
        "Throat mass-flow residual");

    require(
        std::abs(result.continuity_residual) <
            1e-8 * flow_scale,
        "Inlet continuity residual");

    require(
        std::abs(coupled.mass_flow_residual) <
            1e-8 * flow_scale,
        "Coupled mass-flow root residual");

    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    const auto throat = model.evaluate_Tp(
        critical.temperature,
        critical.pressure);

    require(
        static_cast<bool>(inlet) &&
        static_cast<bool>(throat),
        "EOS evaluation failed");

    require(
        std::abs(throat.value.s - inlet.value.s) <
            1e-8,
        "Isentropic entropy conservation");

    // Reconstruct rothalpy directly from the EOS state.
    const double independent_I =
        throat.value.h +
        0.5 * W2 -
        0.5 * result.throat_blade_speed *
            result.throat_blade_speed;

    require(
        std::abs(
            independent_I -
            coupled.inlet_rothalpy) <
            1e-8 * energy_scale,
        "Independent EOS rothalpy check");

    std::cout
        << std::setprecision(12)
        << name
        << " rpm=" << rpm
        << " Wtau=" << Wtau
        << " m/s"
        << " mdot*=" << coupled.critical_mass_flow
        << " kg/s"
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

void test_zero_incidence_equivalence()
{
    thermomodel::IdealGas model;

    const GeometryCoupledChokingInputs aligned_input{
        300.0,
        500000.0,
        9000.0,
        20.0,
        throat_geometry
    };

    const IncidenceCoupledChokingInputs incidence_input{
        300.0,
        500000.0,
        9000.0,
        20.0,
        throat_geometry,
        0.0
    };

    const auto aligned =
        evaluate_geometry_coupled_choking(
            model,
            model,
            impeller_geometry,
            aligned_input);

    const auto incidence =
        evaluate_incidence_coupled_choking(
            model,
            model,
            impeller_geometry,
            incidence_input);

    close(
        incidence.coupled.critical_mass_flow,
        aligned.coupled.critical_mass_flow,
        1e-12,
        "Stage 6C zero-incidence mass-flow equivalence");

    close(
        incidence.throat_relative_velocity.meridional,
        aligned.throat_relative_velocity.meridional,
        1e-12,
        "Stage 6C meridional velocity equivalence");

    close(
        incidence.throat_relative_velocity.circumferential,
        aligned.throat_relative_velocity.circumferential,
        1e-12,
        "Stage 6C circumferential velocity equivalence");

    std::cout
        << "PASS: Stage 6C zero-incidence equivalence\n";
}

void test_sign_symmetry()
{
    thermomodel::IdealGas model;

    const IncidenceCoupledChokingInputs positive{
        300.0,
        500000.0,
        9000.0,
        20.0,
        throat_geometry,
        35.0
    };

    auto negative = positive;
    negative.throat_tangential_relative_velocity = -35.0;

    const auto plus =
        evaluate_incidence_coupled_choking(
            model, model, impeller_geometry, positive);

    const auto minus =
        evaluate_incidence_coupled_choking(
            model, model, impeller_geometry, negative);

    // Under the prescribed-Wtau closure, the
    // thermodynamic equations depend on Wtau^2.
    close(
        plus.coupled.critical_mass_flow,
        minus.coupled.critical_mass_flow,
        1e-12,
        "Prescribed tangential-speed sign symmetry");

    close(
        plus.coupled.throat.critical.sound_speed,
        minus.coupled.throat.critical.sound_speed,
        1e-12,
        "Sound-speed sign symmetry");

    // The coordinate velocities need not be equal,
    // even though the critical thermodynamic states are.
    require(
        std::abs(
            plus.throat_relative_velocity.meridional -
            minus.throat_relative_velocity.meridional) >
            1.0,
        "Tangential sign must affect meridional component");

    std::cout
        << "PASS: prescribed-Wtau sign symmetry\n";
}

void test_invalid_input()
{
    thermomodel::IdealGas model;

    IncidenceCoupledChokingInputs input{
        300.0,
        500000.0,
        9000.0,
        20.0,
        throat_geometry,
        0.0
    };

    input.throat_tangential_relative_velocity =
        std::numeric_limits<double>::quiet_NaN();

    bool rejected = false;

    try {
        (void)evaluate_incidence_coupled_choking(
            model, model, impeller_geometry, input);
    }
    catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Nonfinite Wtau must be rejected");

    std::cout << "PASS: invalid input rejection\n";
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

        test_zero_incidence_equivalence();
        test_sign_symmetry();
        test_invalid_input();

        for (double rpm : {0.0, 9000.0, 15000.0}) {
            for (double Wtau : {-35.0, 0.0, 35.0}) {
                run_case(air, "AirLikeIdeal", rpm, Wtau);
                run_case(
                    methane_ideal, "MethaneIdeal", rpm, Wtau);
                run_case(
                    methane_pr, "MethanePR", rpm, Wtau);
            }
        }

        std::cout
            << "PASS: incidence-coupled choking tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
