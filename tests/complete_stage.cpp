#include <turbocompgpu/complete_stage.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void check(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-9)
{
    const double scale = std::max(1.0, std::abs(expected));

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        std::cerr << "FAIL: " << name
                  << " actual=" << actual
                  << " expected=" << expected << '\n';
        throw std::runtime_error(name);
    }
}

template<class Model>
void test_stationary_reference(
    const Model& model,
    const char* label)
{
    using namespace turbocompgpu;

    const auto inlet_static = evaluate(model, 320.0, 5.0e6);

    const auto inlet = make_stage_station(
        model, inlet_static, 180.0, 2.0);

    const double p_static = inlet_static.p;
    const double p_total =
        inlet.stagnation.stagnation_state.p;

    // Pressure halfway between static and stagnation values.
    const double p_exit =
        p_static + 0.5 * (p_total - p_static);

    const auto result =
        evaluate_lossless_stationary_component(
            model,
            inlet,
            p_exit,
            StationaryComponentKind::Diffuser
        );

    check("mass", result.mass_residual, 0.0);
    check("stagnation enthalpy",
          result.stagnation_enthalpy_residual, 0.0,
          1.0e-5);
    check("entropy", result.entropy_residual, 0.0,
          1.0e-5);
    check("stagnation pressure",
          result.stagnation_pressure_residual, 0.0,
          1.0e-5);

    if (!(result.outlet.absolute_speed <
          result.inlet.absolute_speed)) {
        throw std::runtime_error(
            "Diffuser did not decelerate flow");
    }

    if (!(result.outlet.static_state.p >
          result.inlet.static_state.p)) {
        throw std::runtime_error(
            "Diffuser did not recover static pressure");
    }

    // Check the area returned by continuity.
    check(
        "continuity",
        result.outlet.static_state.rho *
            result.required_effective_area *
            result.outlet.absolute_speed,
        inlet.mass_flow
    );

    // Two stationary components in series.
    const double p_collector =
        p_exit + 0.5 * (p_total - p_exit);

    const auto stage = evaluate_complete_stage_reference(
        model,
        inlet,
        inlet,
        p_exit,
        p_collector
    );

    check("stage mass", stage.mass_residual, 0.0);
    check("stage energy", stage.energy_residual, 0.0,
          1.0e-5);
    check("stage entropy",
          stage.stationary_entropy_change, 0.0,
          1.0e-5);
    check("stage pressure ratio",
          stage.total_to_total_pressure_ratio, 1.0,
          1.0e-5);

    std::cout << "PASS: " << label
              << " lossless stationary stage reference\n";
}

void test_ideal_gas_analytic()
{
    using namespace turbocompgpu;

    thermomodel::IdealGas model;
    const auto& gas = model.parameters();

    const auto state = evaluate(model, 300.0, 5.0e6);
    const auto inlet = make_stage_station(
        model, state, 200.0, 2.0);

    const double T0 =
        300.0 + 200.0 * 200.0 / (2.0 * gas.cp);

    const double p0 =
        5.0e6 * std::pow(T0 / 300.0, gas.cp / gas.R);

    check("analytic T0",
          inlet.stagnation.stagnation_state.T, T0);
    check("analytic p0",
          inlet.stagnation.stagnation_state.p, p0);

    const double p2 = 5.0e6 + 0.4 * (p0 - 5.0e6);

    const auto diffuser =
        evaluate_lossless_stationary_component(
            model,
            inlet,
            p2,
            StationaryComponentKind::Diffuser);

    const double T2 =
        300.0 * std::pow(p2 / 5.0e6, gas.R / gas.cp);

    const double V2 =
        std::sqrt(2.0 * gas.cp * (T0 - T2));

    check("analytic outlet temperature",
          diffuser.outlet.static_state.T, T2);
    check("analytic outlet speed",
          diffuser.outlet.absolute_speed, V2,
          1.0e-7);

    std::cout << "PASS: ideal-gas analytical reference\n";
}

} // namespace

int main()
{
    try {
        test_ideal_gas_analytic();

        thermomodel::IdealGas ideal;
        test_stationary_reference(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters params;
        params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(caloric, params);
        test_stationary_reference(pr, "PengRobinson");

        std::cout << "PASS: Stage 8B.1\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
