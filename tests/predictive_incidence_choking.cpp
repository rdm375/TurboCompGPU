#include <turbocompgpu/predictive_incidence_choking.hpp>

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

using namespace turbocompgpu;

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
            tolerance * std::max(1.0, std::abs(expected)),
        message);
}

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const ThroatGeometry throat{
    12, 0.15, 0.0035, 45.0, 0.002, 0.10
};

template<class Model>
void run_case(
    const Model& model,
    const char* name,
    double rpm)
{
    const PredictiveIncidenceChokingInputs input{
        300.0, 500000.0, rpm, 20.0, throat
    };

    const auto result =
        evaluate_predictive_incidence_choking(
            model, model, geometry, input);

    const auto& choking = result.choking;
    const auto& coupled = choking.coupled;

    const double mdot =
        coupled.critical_mass_flow;

    close(
        result.incidence.inlet_meridional_velocity,
        coupled.inlet_meridional_velocity,
        1e-9,
        "Inlet velocity consistency");

    close(
        result.incidence.inlet_blade_speed,
        coupled.inlet_blade_speed,
        1e-12,
        "Inlet blade-speed consistency");

    close(
        choking.throat_local_velocity.tangential,
        result.incidence.transported_local_velocity.tangential,
        1e-12,
        "Incidence closure");

    require(
        std::abs(result.mass_flow_residual) <
            1e-8 * std::max(1.0, mdot),
        "Outer mass-flow residual");

    require(
        std::abs(result.incidence_residual) < 1e-10,
        "Incidence residual");

    require(
        std::abs(choking.mass_flow_residual) <
            1e-8 * std::max(1.0, mdot),
        "Throat mass-flow residual");

    require(
        std::abs(choking.continuity_residual) <
            1e-8 * std::max(1.0, mdot),
        "Inlet continuity residual");

    require(
        std::abs(coupled.mass_flow_residual) <
            1e-8 * std::max(1.0, mdot),
        "Inner choking root residual");

    require(
        std::abs(choking.sonic_residual) < 1e-10,
        "Sonic normal-velocity residual");

    require(
        std::abs(choking.rothalpy_residual) <
            1e-8 *
            std::max(1.0, std::abs(coupled.inlet_rothalpy)),
        "Rothalpy conservation");

    require(
        choking.throat_relative_velocity.meridional > 0.0,
        "Nonpositive meridional throat velocity");

    // Independent direct EOS check.
    const auto inlet = model.evaluate_Tp(
        input.inlet_temperature,
        input.inlet_pressure);

    const auto critical = model.evaluate_Tp(
        coupled.throat.critical.temperature,
        coupled.throat.critical.pressure);

    require(
        static_cast<bool>(inlet) &&
        static_cast<bool>(critical),
        "Independent EOS evaluation");

    require(
        std::abs(critical.value.s - inlet.value.s) < 1e-8,
        "Independent entropy conservation");

    const auto& velocity =
        choking.throat_relative_velocity;

    const double independent_rothalpy =
        critical.value.h +
        0.5 * (
            velocity.meridional * velocity.meridional +
            velocity.circumferential * velocity.circumferential
        ) -
        0.5 * choking.throat_blade_speed *
            choking.throat_blade_speed;

    require(
        std::abs(
            independent_rothalpy -
            coupled.inlet_rothalpy) <
            1e-8 *
            std::max(1.0, std::abs(coupled.inlet_rothalpy)),
        "Independent rothalpy conservation");

    // Independently re-evaluate the Stage 6D capacity
    // using the converged projected tangential velocity.
    const auto independently_coupled =
        evaluate_incidence_coupled_choking(
            model,
            model,
            geometry,
            {
                input.inlet_temperature,
                input.inlet_pressure,
                input.rotational_speed_rpm,
                input.inlet_prewhirl,
                input.throat_geometry,
                result.incidence
                    .transported_local_velocity.tangential
            });

    close(
        independently_coupled.coupled.critical_mass_flow,
        mdot,
        1e-10,
        "Stage 6D capacity consistency");

    std::cout
        << std::setprecision(12)
        << name
        << " rpm=" << rpm
        << " mdot*=" << mdot
        << " kg/s"
        << " Wtau="
        << choking.throat_local_velocity.tangential
        << " m/s"
        << " incidence="
        << result.incidence.incidence_angle_deg
        << " deg"
        << " F=" << result.mass_flow_residual
        << " kg/s"
        << " dI=" << choking.rothalpy_residual
        << " J/kg"
        << " iterations=" << result.iterations
        << '\n';
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

        for (double rpm : {0.0, 9000.0, 12000.0, 15000.0}) {
            run_case(air, "AirLikeIdeal", rpm);
            run_case(methane_ideal, "MethaneIdeal", rpm);
            run_case(methane_pr, "MethanePR", rpm);
        }

        std::cout
            << "PASS: predictive-incidence choking tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
