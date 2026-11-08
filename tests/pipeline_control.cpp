#include <turbocompgpu/pipeline_compressor.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/caloric.hpp>
#include <thermomodel/peng_robinson.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numbers>
#include <memory>
#include <limits>
#include <utility>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const char* label)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance * scale) {
        std::cerr << label << ": actual=" << actual
                  << " expected=" << expected << '\n';
        throw std::runtime_error(label);
    }
}

template<class Exception, class Function>
void throws(Function&& function, const char* label)
{
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(label);
}

struct ComponentFixture {
    thermomodel::IdealGas model;
    ImpellerSystemGeometry impeller{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
    VanelessDiffuserGeometry diffuser;
    CompressorComponentConfig config;

    ComponentFixture()
    {
        const double r2 = impeller.exit.diameter / 2.0;
        const double area = effective_flow_area(
            impeller_exit_area(impeller.exit), impeller.exit.blockage);
        diffuser = {r2, 1.5 * r2, area / (2.0 * std::numbers::pi * r2), 0.0};
        config.slip_factor = 0.90;
        config.residual_scales = {50000.0, 100.0, 500000.0};
        config.solver_options.residual_tolerance = 1.0e-13;
        config.initial_outlet_temperature_k = 340.0;
        config.initial_outlet_pressure_pa = 700000.0;
    }

    CompressorComponent make(const GuideVaneMap* map = nullptr) const
    {
        return CompressorComponent(
            model, model, impeller, diffuser, config, map);
    }
};

void test_component_contract_and_vanes()
{
    ComponentFixture f;
    auto fixed = f.make();
    CompressorComponentInput input{
        {300.0, 500000.0}, 5.0, 12000.0, std::nullopt};
    const auto baseline = fixed.evaluate(input);
    require(baseline.valid, "fixed-vane component solve failed");
    require(baseline.outlet_static_state.has_value(), "missing outlet static state");
    require(baseline.outlet_total_state.has_value(), "missing outlet total state");
    require(baseline.shaft_power_w > 0.0, "nonpositive shaft power");
    close(baseline.compressor_throughflow_kg_s, 5.0, 1.0e-12, "through-flow");
    close(baseline.stage.stage->mass_flow, 5.0, 1.0e-12, "stage mass flow");
    close(baseline.stage.stage->mass_residual, 0.0, 1.0e-10,
          "component mass residual");
    const auto& station3 = *baseline.stage.stage->station3;
    const double outlet_kinetic_energy = 0.5 * (
        station3.meridional_velocity * station3.meridional_velocity +
        station3.tangential_velocity * station3.tangential_velocity);
    close(baseline.outlet_total_state->h - baseline.outlet_static_state->h,
          outlet_kinetic_energy, 1.0e-10,
          "static-to-total outlet enthalpy");
    require(baseline.outlet_total_state->p >= baseline.outlet_static_state->p,
            "total outlet pressure below static pressure");
    close(baseline.shaft_power_w,
          baseline.stage.stage->impeller_shaft_power, 1.0e-12,
          "component shaft power");

    GuideVaneMap map;
    map.domain = {0.0, 30.0, 500.0, 20000.0, 0.0, 20.0};
    map.calibrated_effects = [](double, double, double,
                                const thermomodel::State&) {
        return GuideVaneEffects{0.0, 0.0};
    };
    const auto same_as_fixed = f.make(&map);
    input.guide_vane_angle_deg = 15.0;
    const auto mapped = same_as_fixed.evaluate(input);
    require(mapped.valid, "mapped fixed-vane case failed");
    close(mapped.outlet_static_state->p,
          baseline.outlet_static_state->p, 1.0e-9,
          "zero-effect vane map pressure");
    close(mapped.shaft_power_w, baseline.shaft_power_w, 1.0e-9,
          "zero-effect vane map power");

    throws<std::out_of_range>([&] {
        (void)map.evaluate(31.0, 12000.0, 5.0, {});
    }, "out-of-domain guide-vane query accepted");
    map.calibrated_effects = [](double, double, double,
                                const thermomodel::State&) {
        return GuideVaneEffects{0.0, -1.0};
    };
    throws<std::runtime_error>([&] {
        (void)map.evaluate(15.0, 12000.0, 5.0, {});
    }, "invalid guide-vane calibration accepted");

    input.guide_vane_angle_deg = 15.0;
    throws<std::invalid_argument>([&] { (void)fixed.evaluate(input); },
                                  "guide-vane command without map accepted");
    input.guide_vane_angle_deg.reset();
    input.compressor_throughflow_kg_s = 0.0;
    throws<std::invalid_argument>([&] { (void)fixed.evaluate(input); },
                                  "invalid component flow accepted");

    auto bad_diffuser = f.diffuser;
    bad_diffuser.inlet_radius *= 1.1;
    CompressorComponent invalid_component(
        f.model, f.model, f.impeller, bad_diffuser, f.config);
    input.compressor_throughflow_kg_s = 5.0;
    const auto failed = invalid_component.evaluate(input);
    require(!failed.valid && !failed.message.empty(),
            "diffuser failure was not returned through the component contract");
}

SteadyPlantMeasurements plant(double rpm)
{
    return {
        350000.0 + 50.0 * rpm,
        500000.0 + 200.0 * rpm,
        0.003 * rpm,
        0.005 * rpm * rpm
    };
}

void test_steady_control_modes_and_limits()
{
    SteadySetpointOptions limits;
    limits.minimum_speed_rpm = 1000.0;
    limits.maximum_speed_rpm = 2000.0;
    limits.maximum_shaft_power_w = 11250.0;
    limits.absolute_tolerance = 1.0e-5;
    limits.relative_tolerance = 1.0e-9;

    for (const auto [mode, setpoint] : {
             std::pair{SteadySetpointMode::SuctionPressure, 425000.0},
             std::pair{SteadySetpointMode::DischargePressure, 800000.0},
             std::pair{SteadySetpointMode::PipelineExportFlow, 4.5}}) {
        const auto result = solve_steady_speed_setpoint(
            mode, setpoint, limits, plant);
        require(result.status == SteadySetpointStatus::Converged,
                "automatic setpoint mode failed");
        close(result.speed_rpm, 1500.0, 1.0e-6,
              "automatic setpoint speed");
        require(result.measurements.compressor_shaft_power_w <=
                    limits.maximum_shaft_power_w + 1.0e-6,
                "automatic control exceeded shaft-power limit");
    }

    const auto unreachable = solve_steady_speed_setpoint(
        SteadySetpointMode::DischargePressure, 900000.0, limits, plant);
    require(unreachable.status == SteadySetpointStatus::PowerLimitReached,
            "unreachable pressure did not report power saturation");
    close(unreachable.speed_rpm, 1500.0, 1.0e-6,
          "power-limited speed");

    auto manual_limits = limits;
    manual_limits.maximum_shaft_power_w =
        std::numeric_limits<double>::infinity();
    const auto manual = solve_steady_speed_setpoint(
        SteadySetpointMode::ManualSpeed, 2500.0, manual_limits, plant);
    require(manual.status == SteadySetpointStatus::SpeedLimitReached,
            "manual speed clamp not reported");
    close(manual.speed_rpm, 2000.0, 1.0e-12, "manual speed clamp");

    auto unbounded = limits;
    unbounded.maximum_shaft_power_w =
        std::numeric_limits<double>::infinity();
    const auto not_bracketed = solve_steady_speed_setpoint(
        SteadySetpointMode::DischargePressure, 950000.0, unbounded, plant);
    require(not_bracketed.status == SteadySetpointStatus::TargetNotBracketed,
            "unbracketed setpoint status incorrect");

    auto one_iteration = unbounded;
    one_iteration.maximum_iterations = 1;
    const auto limited_iterations = solve_steady_speed_setpoint(
        SteadySetpointMode::DischargePressure, 775000.0,
        one_iteration, plant);
    require(limited_iterations.status == SteadySetpointStatus::IterationLimit,
            "root iteration limit status incorrect");

    auto infeasible_power = limits;
    infeasible_power.maximum_shaft_power_w = 4000.0;
    const auto no_feasible_speed = solve_steady_speed_setpoint(
        SteadySetpointMode::DischargePressure, 800000.0,
        infeasible_power, plant);
    require(no_feasible_speed.status == SteadySetpointStatus::PowerLimitReached,
            "power limit below minimum-speed demand was not reported");

    const auto manual_power_limited = solve_steady_speed_setpoint(
        SteadySetpointMode::ManualSpeed, 2000.0, limits, plant);
    require(manual_power_limited.status ==
                SteadySetpointStatus::PowerLimitReached,
            "manual power limit not reported");
    close(manual_power_limited.speed_rpm, 1500.0, 1.0e-6,
          "manual power-limited speed");

    throws<std::runtime_error>([&] {
        (void)solve_steady_speed_setpoint(
            SteadySetpointMode::DischargePressure, 800000.0, limits,
            [](double) {
                return SteadyPlantMeasurements{1.0, 2.0, 3.0,
                                               std::numeric_limits<double>::quiet_NaN()};
            });
    }, "invalid plant measurements accepted");
}

void test_recycle_thermodynamics_and_loop()
{
    thermomodel::IdealGas model;
    const auto fresh = evaluate(model, 300.0, 500000.0);
    const auto hot = evaluate(model, 330.0, 500000.0);
    const CompressorPipePort fresh_port{300.0, 500000.0};

    const auto no_recycle = mix_suction_streams(
        model, fresh_port, 5.0, hot, 0.0);
    close(no_recycle.T, fresh.T, 1.0e-12, "zero-recycle temperature");
    close(no_recycle.h, fresh.h, 1.0e-12, "zero-recycle enthalpy");

    const auto mixed = mix_suction_streams(
        model, fresh_port, 5.0, hot, 2.0);
    const double expected_h = (5.0 * fresh.h + 2.0 * hot.h) / 7.0;
    close(mixed.h, expected_h, 1.0e-10, "suction enthalpy balance");
    close(mixed.p, fresh_port.static_pressure_pa, 1.0e-12,
          "suction mixer pressure");

    const auto discharge = evaluate(model, 350.0, 1000000.0);
    const auto returned = isenthalpic_recycle_return(
        model, discharge, 500000.0);
    close(returned.h, discharge.h, 1.0e-10,
          "isenthalpic recycle return");
    close(returned.p, 500000.0, 1.0e-12,
          "recycle return pressure");

    RecycleValveMap valve;
    valve.domain = {400000.0, 600000.0, 900000.0, 1100000.0, 2.0};
    valve.calibrated_flow_kg_s = [](double opening, double suction,
                                   const thermomodel::State& outlet) {
        return opening * (outlet.p - suction) / 1000000.0;
    };
    close(valve.flow(0.0, 500000.0, discharge), 0.0, 1.0e-12,
          "closed recycle valve");
    close(valve.flow(0.5, 500000.0, discharge), 0.25, 1.0e-12,
          "recycle pressure-drop flow law");
    throws<std::out_of_range>([&] {
        (void)valve.flow(0.5, 700000.0, discharge);
    }, "out-of-domain recycle query accepted");
    auto excessive_valve = valve;
    excessive_valve.calibrated_flow_kg_s = [](double, double,
                                               const thermomodel::State&) {
        return 3.0;
    };
    throws<std::runtime_error>([&] {
        (void)excessive_valve.flow(0.5, 500000.0, discharge);
    }, "recycle valve flow above its declared maximum accepted");

    auto compressor = [&model](const thermomodel::State& suction,
                               double throughflow) {
        CompressorComponentResult result;
        result.valid = true;
        result.compressor_throughflow_kg_s = throughflow;
        result.outlet_static_state = evaluate(
            model, suction.T + 30.0, 1000000.0);
        result.shaft_power_w = 100.0 * throughflow;
        return result;
    };
    auto return_model = [&model](const thermomodel::State& outlet,
                                 double suction_pressure) {
        return isenthalpic_recycle_return(
            model, outlet, suction_pressure);
    };
    RecycleSteadyOptions options;
    options.maximum_iterations = 80;
    options.relaxation = 1.0;
    options.relative_return_enthalpy_tolerance = 1.0e-6;

    const auto coupled = solve_steady_recycle(
        model, fresh_port, 5.0, 0.5, valve,
        compressor, return_model, options);
    require(coupled.converged, "steady recycle iteration did not converge");
    close(coupled.recycle_flow_kg_s, 0.25, 1.0e-8,
          "recycle valve flow");
    close(coupled.compressor_throughflow_kg_s, 5.25, 1.0e-8,
          "compressor through-flow accounting");
    close(coupled.recycle_mass_residual_kg_s, 0.0, 1.0e-8,
          "recycle mass residual");
    close(coupled.suction_mixer_enthalpy_residual_w, 0.0, 1.0e-6,
          "suction mixer energy residual");
    close(coupled.recycle_cooler_heat_removed_w, 0.0, 1.0e-6,
          "uncooled recycle duty");

    auto cooler_return = [&model](const thermomodel::State& outlet,
                                  double suction_pressure) {
        const auto throttled = isenthalpic_recycle_return(
            model, outlet, suction_pressure);
        return evaluate(model, throttled.T - 10.0, suction_pressure);
    };
    const auto cooled = solve_steady_recycle(
        model, fresh_port, 5.0, 0.5, valve,
        compressor, cooler_return, options);
    require(cooled.converged, "cooled recycle iteration did not converge");
    require(cooled.recycle_cooler_heat_removed_w > 0.0,
            "aftercooler duty was not reported positive");

    auto invalid_return = [&model](const thermomodel::State& outlet,
                                   double suction_pressure) {
        return evaluate(model, outlet.T, suction_pressure * 0.9);
    };
    throws<std::runtime_error>([&] {
        (void)solve_steady_recycle(
            model, fresh_port, 5.0, 0.5, valve,
            compressor, invalid_return, options);
    }, "recycle return at the wrong pressure was accepted");

    auto one_step = options;
    one_step.maximum_iterations = 1;
    const auto nonconverged = solve_steady_recycle(
        model, fresh_port, 5.0, 0.5, valve,
        compressor, return_model, one_step);
    require(!nonconverged.converged,
            "recycle iteration limit did not report failure");

    const auto closed = solve_steady_recycle(
        model, fresh_port, 5.0, 0.0, valve,
        compressor, return_model, options);
    require(closed.converged, "closed recycle case did not converge");
    close(closed.recycle_flow_kg_s, 0.0, 1.0e-12,
          "closed recycle flow");
    close(closed.compressor_throughflow_kg_s, 5.0, 1.0e-12,
          "closed-valve compressor flow");

    auto caloric = std::make_shared<thermomodel::ConstantCp>(2200.0);
    thermomodel::PRParameters parameters;
    parameters.root_policy = thermomodel::PRRootPolicy::Vapor;
    thermomodel::PengRobinson peng_robinson(caloric, parameters);
    const auto pr_fresh = evaluate(peng_robinson, 300.0, 500000.0);
    const auto pr_return = evaluate(peng_robinson, 330.0, 500000.0);
    const auto pr_mix = mix_suction_streams(
        peng_robinson, fresh_port, 5.0, pr_return, 2.0);
    close(pr_mix.h,
          (5.0 * pr_fresh.h + 2.0 * pr_return.h) / 7.0,
          1.0e-9, "Peng-Robinson suction enthalpy balance");
    close(pr_mix.p, fresh_port.static_pressure_pa, 1.0e-12,
          "Peng-Robinson suction pressure");
}

} // namespace

int main()
{
    try {
        test_component_contract_and_vanes();
        test_steady_control_modes_and_limits();
        test_recycle_thermodynamics_and_loop();
        std::cout << "Pipeline compressor/control/recycle tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
