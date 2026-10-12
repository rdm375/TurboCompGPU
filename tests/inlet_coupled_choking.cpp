#include <turbocompgpu/inlet_coupled_choking.hpp>

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

template<class Model>
void run_case(
    const Model& model,
    const char* name,
    double rpm)
{
    const turbocompgpu::ImpellerSystemGeometry geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const turbocompgpu::InletCoupledChokingInputs input{
        300.0,
        500000.0,
        rpm,
        20.0,
        0.15,
        -35.0,
        0.002
    };

    const auto result =
        turbocompgpu::evaluate_inlet_coupled_choking(
            model, model, geometry, input);

    require(
        result.critical_mass_flow > 0.0,
        "Nonpositive critical mass flow");

    require(
        std::abs(result.mass_flow_residual) <
            1e-8 * std::max(1.0, result.critical_mass_flow),
        "Mass-flow coupling residual too large");

    const auto inlet =
        model.evaluate_Tp(300.0, 500000.0);

    require(static_cast<bool>(inlet), "Inlet EOS failed");

    const auto effective = model.evaluate_Tp(
        result.effective_stagnation_temperature,
        result.effective_stagnation_pressure);

    require(static_cast<bool>(effective), "Effective EOS failed");

    const double expected_h =
        result.inlet_rothalpy +
        0.5 * (
            result.throat.blade_speed *
                result.throat.blade_speed -
            input.throat_tangential_relative_velocity *
                input.throat_tangential_relative_velocity);

    const double energy_scale =
        std::max(1.0, std::abs(expected_h));

    require(
        std::abs(effective.value.h - expected_h) <
            1e-8 * energy_scale,
        "Effective stagnation enthalpy mismatch");

    require(
        std::abs(effective.value.s - inlet.value.s) < 1e-8,
        "Effective stagnation entropy mismatch");

    std::cout
        << std::setprecision(12)
        << name
        << " rpm=" << rpm
        << " mdot*=" << result.critical_mass_flow
        << " kg/s"
        << " T0eff=" << result.effective_stagnation_temperature
        << " K"
        << " p0eff=" << result.effective_stagnation_pressure
        << " Pa"
        << " residual=" << result.mass_flow_residual
        << " kg/s\n";
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

        for (double rpm : {0.0, 9000.0, 12000.0, 15000.0}) {
            run_case(air, "AirLikeIdeal", rpm);
            run_case(methane_ideal, "MethaneIdeal", rpm);
            run_case(methane_pr, "MethanePR", rpm);
        }

        std::cout
            << "PASS: inlet-coupled choking tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
