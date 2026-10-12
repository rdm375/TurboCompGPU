#include <turbocompgpu/passage_turning_choking.hpp>
#include <turbocompgpu/geometry_coupled_choking.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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
    const auto aligned =
        evaluate_geometry_coupled_choking(
            model, model, geometry,
            {300.0, 500000.0, rpm, 20.0, throat});

    const auto frozen =
        evaluate_predictive_incidence_choking(
            model, model, geometry,
            {300.0, 500000.0, rpm, 20.0, throat});

    double previous_mdot =
        std::numeric_limits<double>::quiet_NaN();

    for (double lambda : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto result =
            evaluate_passage_turning_choking(
                model, model, geometry,
                {
                    300.0,
                    500000.0,
                    rpm,
                    20.0,
                    throat,
                    lambda
                });

        const auto& choking = result.choking;
        const double mdot =
            choking.coupled.critical_mass_flow;

        close(
            choking.throat_local_velocity.tangential,
            lambda *
                result.incidence
                    .transported_local_velocity.tangential,
            1e-12,
            "Turning closure");

        require(
            std::abs(result.turning_residual) < 1e-10,
            "Turning residual");

        require(
            std::abs(result.mass_flow_residual) <
                1e-8 * std::max(1.0, mdot),
            "Outer mass-flow residual");

        require(
            std::abs(choking.rothalpy_residual) <
                1e-8 *
                std::max(
                    1.0,
                    std::abs(choking.coupled.inlet_rothalpy)),
            "Rothalpy residual");

        require(
            std::abs(choking.continuity_residual) <
                1e-8 * std::max(1.0, mdot),
            "Continuity residual");

        require(
            std::abs(choking.sonic_residual) < 1e-10,
            "Sonic residual");

        require(
            choking.throat_relative_velocity.meridional > 0.0,
            "Forward meridional flow");

        if (lambda == 0.0) {
            close(
                mdot,
                aligned.coupled.critical_mass_flow,
                1e-9,
                "Stage 6C endpoint mass flow");

            close(
                choking.throat_local_velocity.tangential,
                0.0,
                1e-12,
                "Stage 6C endpoint tangential velocity");
        }

        if (lambda == 1.0) {
            close(
                mdot,
                frozen.choking.coupled.critical_mass_flow,
                1e-9,
                "Stage 6E endpoint mass flow");

            close(
                choking.throat_local_velocity.tangential,
                frozen.choking.throat_local_velocity.tangential,
                1e-9,
                "Stage 6E endpoint tangential velocity");
        }

        // With this closure, |Wtau| can change nonmonotonically
        // as lambda varies because mdot is solved simultaneously.
        // Do not impose monotonic mass-flow behavior.
        previous_mdot = mdot;

        std::cout
            << std::setprecision(12)
            << name
            << " rpm=" << rpm
            << " lambda=" << lambda
            << " mdot*=" << mdot
            << " kg/s Wtau="
            << choking.throat_local_velocity.tangential
            << " m/s"
            << " F=" << result.mass_flow_residual
            << " kg/s"
            << " dI=" << choking.rothalpy_residual
            << " J/kg"
            << " iterations=" << result.iterations
            << '\n';
    }

    (void)previous_mdot;
}

void test_invalid_lambda()
{
    thermomodel::IdealGas air;

    for (double lambda : {
             -0.01,
             1.01,
             std::numeric_limits<double>::quiet_NaN()
         }) {
        bool rejected = false;

        try {
            (void)evaluate_passage_turning_choking(
                air, air, geometry,
                {300.0, 500000.0, 9000.0, 20.0, throat, lambda});
        }
        catch (const std::invalid_argument&) {
            rejected = true;
        }

        require(rejected, "Invalid lambda not rejected");
    }

    std::cout << "PASS: invalid turning-fraction rejection\n";
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

        thermomodel::IdealGas methane_ideal(methane_params);

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PengRobinson methane_pr(
            caloric, pr_params);

        for (double rpm : {0.0, 9000.0, 15000.0}) {
            run_case(air, "AirLikeIdeal", rpm);
            run_case(methane_ideal, "MethaneIdeal", rpm);
            run_case(methane_pr, "MethanePR", rpm);
        }

        test_invalid_lambda();

        std::cout
            << "PASS: passage-turning choking tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
