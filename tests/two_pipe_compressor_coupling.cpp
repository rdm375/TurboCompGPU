#include <turbocompgpu/transient_pipe_coupling.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <pipe_sim/physics.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

struct Fixture {
    thermomodel::IdealGas model;
    ImpellerSystemGeometry impeller{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
    VanelessDiffuserGeometry diffuser;
    CompressorComponentConfig compressor_config;

    Fixture()
    {
        const double r2 = impeller.exit.diameter / 2.0;
        const double area = effective_flow_area(
            impeller_exit_area(impeller.exit), impeller.exit.blockage);
        diffuser = {
            r2, 1.5 * r2,
            area / (2.0 * std::numbers::pi * r2), 0.0
        };
        compressor_config.slip_factor = 0.90;
        compressor_config.residual_scales = {50000.0, 100.0, 500000.0};
        compressor_config.solver_options.residual_tolerance = 1.0e-13;
        compressor_config.initial_outlet_temperature_k = 340.0;
        compressor_config.initial_outlet_pressure_pa = 700000.0;
    }

    CompressorComponent compressor() const
    {
        return CompressorComponent(
            model, model, impeller, diffuser, compressor_config);
    }
};

PipeSegmentModel segment_model()
{
    constexpr double diameter = 0.60;
    const double area = std::numbers::pi * diameter * diameter / 4.0;
    return {
        2,
        {1000.0, 1.0, 0.65, diameter, area, 300.0,
         1.0, 287.0, 1.8e-5, 1.0e-5},
        {1.0e-8, 1.0e-12, 40, false}
    };
}

PipeSegmentState steady_segment_state(
    const PipeSegmentModel& model,
    double inlet_pressure_pa,
    double flow_kg_s)
{
    const auto& par = model.parameters;
    const double reynolds = pipe_sim::reynolds_mass(
        flow_kg_s, par.diameter, par.area, par.viscosity);
    const double friction = pipe_sim::friction_sj(
        reynolds, par.roughness, par.diameter);
    const double coefficient =
        friction * par.z * par.gas_constant * par.temperature *
        flow_kg_s * flow_kg_s /
        (par.diameter * par.area * par.area);

    PipeSegmentState state;
    state.pressure_pa.resize(static_cast<std::size_t>(model.cells + 1));
    state.flow_kg_s.assign(static_cast<std::size_t>(model.cells), flow_kg_s);
    state.inlet_flow_kg_s = flow_kg_s;
    state.outlet_flow_kg_s = flow_kg_s;
    for (int i = 0; i <= model.cells; ++i) {
        const double pressure_squared = inlet_pressure_pa * inlet_pressure_pa -
            coefficient * static_cast<double>(i) * par.dx;
        if (!(pressure_squared > 0.0)) {
            throw std::runtime_error("Invalid steady pipe initializer");
        }
        state.pressure_pa[static_cast<std::size_t>(i)] =
            std::sqrt(pressure_squared);
    }
    return state;
}

double pipe_linepack(
    const PipeSegmentModel& model,
    const PipeSegmentState& state)
{
    double density_sum = 0.0;
    for (int i = 0; i <= model.cells; ++i) {
        const double weight =
            (i == 0 || i == model.cells) ? 0.5 : 1.0;
        density_sum += weight * pipe_sim::density_cz(
            state.pressure_pa[static_cast<std::size_t>(i)],
            model.parameters.temperature,
            model.parameters.z,
            model.parameters.gas_constant);
    }
    return model.parameters.area * model.parameters.dx * density_sum;
}

double pipe_mass_balance_defect(
    const PipeSegmentModel& model,
    const PipeSegmentState& old_state,
    const PipeSegmentState& new_state)
{
    const double theta = model.parameters.theta;
    const double inventory_rate =
        (pipe_linepack(model, new_state) -
         pipe_linepack(model, old_state)) / model.parameters.dt;
    const double boundary_rate =
        theta * (new_state.inlet_flow_kg_s - new_state.outlet_flow_kg_s) +
        (1.0 - theta) *
            (old_state.inlet_flow_kg_s - old_state.outlet_flow_kg_s);
    return inventory_rate - boundary_rate;
}

void test_two_pipe_compressor_step()
{
    Fixture f;
    const auto compressor = f.compressor();
    const auto pipe = segment_model();
    TwoPipeCompressorStepOptions coupling_options;
    coupling_options.minimum_compressor_throughflow_kg_s = 4.0;
    coupling_options.maximum_compressor_throughflow_kg_s = 6.0;
    const auto upstream_initial = steady_segment_state(
        pipe, 500000.0, 5.0);
    const auto initial_compressor = compressor.evaluate({
        {300.0, upstream_initial.pressure_pa.back()},
        5.0, 12000.0, std::nullopt
    });
    require(initial_compressor.valid &&
                initial_compressor.outlet_static_state.has_value(),
            "initial compressor operating point did not converge");
    const auto downstream_initial = steady_segment_state(
        pipe, initial_compressor.outlet_static_state->p, 5.0);
    const auto direct = advance_two_pipe_compressor(
        pipe, upstream_initial, 500000.0,
        pipe, downstream_initial, 5.0,
        compressor, 12000.0, coupling_options);

    require(direct.valid, direct.message.c_str());

    RecycleValveMap valve_map;
    valve_map.domain = {100000.0, 900000.0, 100000.0, 2000000.0, 1.0};
    valve_map.calibrated_flow_kg_s = [](
        double opening, double suction_pressure,
        const thermomodel::State& discharge) {
        return opening * (discharge.p - suction_pressure) / 1.0e6;
    };
    TwoPipeRecycleCouplingOptions recycle_options;
    recycle_options.minimum_fresh_flow_kg_s = 4.0;
    recycle_options.maximum_fresh_flow_kg_s = 6.0;
    recycle_options.minimum_compressor_throughflow_kg_s = 4.0;
    recycle_options.maximum_compressor_throughflow_kg_s = 7.0;
    recycle_options.recycle_loop.maximum_iterations = 80;
    const auto recycle_return = [&f](
        const thermomodel::State& discharge, double suction_pressure) {
        return isenthalpic_recycle_return(
            f.model, discharge, suction_pressure);
    };
    const auto closed_recycle = advance_two_pipe_compressor_with_recycle(
        f.model, pipe, upstream_initial, 500000.0,
        pipe, downstream_initial, 5.0,
        compressor, 12000.0, 0.0, valve_map, recycle_return,
        recycle_options);
    require(closed_recycle.valid, closed_recycle.message.c_str());
    require(std::abs(closed_recycle.recycle_flow_kg_s) < 1.0e-12,
            "closed recycle valve passed nonzero flow");
    require(std::abs(
                closed_recycle.compressor_throughflow_kg_s -
                direct.compressor.compressor_throughflow_kg_s) < 1.0e-6,
            "closed recycle case did not recover the no-recycle solution");

    const auto recycled = advance_two_pipe_compressor_with_recycle(
        f.model, pipe, upstream_initial, 500000.0,
        pipe, downstream_initial, 5.0,
        compressor, 12000.0, 0.5, valve_map, recycle_return,
        recycle_options);
    require(recycled.valid, recycled.message.c_str());
    require(recycled.recycle_flow_kg_s > 0.0,
            "open recycle valve produced no recycle flow");
    require(recycled.compressor_throughflow_kg_s >
                recycled.fresh_flow_kg_s,
            "compressor through-flow did not include recycle flow");
    require(std::abs(recycled.compressor_node_mass_residual_kg_s) < 1.0e-6,
            "recycle compressor-node mass balance failed");
    require(std::abs(recycled.recycle_mass_residual_kg_s) < 1.0e-6,
            "recycle valve-flow iteration did not close");
    require(std::abs(recycled.suction_mixer_enthalpy_residual_w) < 1.0e-6,
            "recycle suction-mixer energy balance failed");
    require(std::abs(pipe_mass_balance_defect(
                pipe, upstream_initial, recycled.upstream)) < 1.0e-6,
            "recycle upstream linepack balance failed");
    require(std::abs(pipe_mass_balance_defect(
                pipe, downstream_initial, recycled.downstream)) < 1.0e-6,
            "recycle downstream linepack balance failed");

    RecycleDeliveryPressureBoundaryOptions recycle_delivery_options;
    recycle_delivery_options.recycle_coupling = recycle_options;
    recycle_delivery_options.minimum_export_flow_kg_s = 4.0;
    recycle_delivery_options.maximum_export_flow_kg_s = 6.0;
    recycle_delivery_options.absolute_pressure_tolerance_pa = 1.0e-4;
    recycle_delivery_options.relative_pressure_tolerance = 1.0e-10;
    const auto recycled_pressure_boundary =
        advance_two_pipe_compressor_with_recycle_at_delivery_pressure(
            f.model, pipe, upstream_initial, 500000.0,
            pipe, downstream_initial,
            recycled.measurements.discharge_pressure_pa,
            compressor, 12000.0, 0.5, valve_map, recycle_return,
            recycle_delivery_options);
    require(recycled_pressure_boundary.valid,
            recycled_pressure_boundary.message.c_str());
    require(std::abs(recycled_pressure_boundary.measurements.pipeline_export_flow_kg_s -
                     5.0) < 1.0e-6,
            "recycle delivery-pressure solve did not recover the expected net export");
    require(recycled_pressure_boundary.recycle_flow_kg_s > 0.0,
            "delivery-pressure recycle solve lost recycle flow");
    require(std::abs(recycled_pressure_boundary.compressor_node_mass_residual_kg_s) <
                1.0e-6,
            "delivery-pressure recycle compressor-node balance failed");

    SteadySetpointOptions recycled_flow_control_options;
    recycled_flow_control_options.minimum_speed_rpm = 12000.0;
    recycled_flow_control_options.maximum_speed_rpm = 12001.0;
    recycled_flow_control_options.absolute_tolerance = 1.0e-8;
    recycled_flow_control_options.relative_tolerance = 1.0e-10;
    const auto recycle_flow_controlled =
        advance_flow_controlled_two_pipe_compressor_with_recycle_at_delivery_pressure(
            recycled_pressure_boundary.measurements.pipeline_export_flow_kg_s,
            recycled_flow_control_options,
            f.model, pipe, upstream_initial, 500000.0,
            pipe, downstream_initial,
            recycled.measurements.discharge_pressure_pa,
            compressor, 0.5, valve_map, recycle_return,
            recycle_delivery_options);
    require(recycle_flow_controlled.control.status ==
                SteadySetpointStatus::Converged,
            "net-flow controller did not converge with recycle enabled");
    require(recycle_flow_controlled.network.valid,
            recycle_flow_controlled.network.message.c_str());
    require(std::abs(
                recycle_flow_controlled.network.measurements.pipeline_export_flow_kg_s -
                recycled_pressure_boundary.measurements.pipeline_export_flow_kg_s) <
                1.0e-5,
            "recycle flow controller missed the net-export setpoint");
    require(std::abs(
                recycle_flow_controlled.network.compressor_throughflow_kg_s -
                recycle_flow_controlled.network.fresh_flow_kg_s -
                recycle_flow_controlled.network.recycle_flow_kg_s) < 1.0e-8,
            "recycle flow controller did not account for fresh feed plus recycle");
    require(recycle_flow_controlled.network.compressor_throughflow_kg_s >
                recycle_flow_controlled.network.fresh_flow_kg_s,
            "recycle flow controller did not distinguish through-flow from fresh feed");

    DeliveryPressureBoundaryOptions delivery_options;
    delivery_options.compressor_coupling = coupling_options;
    delivery_options.minimum_export_flow_kg_s = 4.0;
    delivery_options.maximum_export_flow_kg_s = 6.0;
    delivery_options.absolute_pressure_tolerance_pa = 1.0e-4;
    delivery_options.relative_pressure_tolerance = 1.0e-10;
    const double delivery_pressure_pa = downstream_initial.pressure_pa.back();
    const auto pressure_boundary_step =
        advance_two_pipe_compressor_at_delivery_pressure(
            pipe, upstream_initial, 500000.0,
            pipe, downstream_initial, delivery_pressure_pa,
            compressor, 12000.0, delivery_options);
    require(pressure_boundary_step.valid,
            pressure_boundary_step.message.c_str());
    require(std::abs(pressure_boundary_step.measurements.discharge_pressure_pa -
                     delivery_pressure_pa) < 2.0,
            "downstream delivery pressure boundary was not satisfied");

    SteadySetpointOptions flow_control_options;
    flow_control_options.minimum_speed_rpm = 11990.0;
    flow_control_options.maximum_speed_rpm = 12010.0;
    flow_control_options.absolute_tolerance = 1.0e-9;
    flow_control_options.relative_tolerance = 1.0e-10;
    flow_control_options.maximum_iterations = 60;
    const auto flow_controlled =
        advance_flow_controlled_two_pipe_compressor_at_delivery_pressure(
            pressure_boundary_step.measurements.pipeline_export_flow_kg_s,
            flow_control_options,
            pipe, upstream_initial, 500000.0,
            pipe, downstream_initial, delivery_pressure_pa,
            compressor, delivery_options);
    require(flow_controlled.control.status == SteadySetpointStatus::Converged,
            "pipeline-flow controller did not converge against calculated export flow");
    require(flow_controlled.network.valid,
            flow_controlled.network.message.c_str());
    require(std::abs(flow_controlled.control.speed_rpm - 12000.0) < 0.5,
            "pipeline-flow controller did not recover the target compressor speed");
    require(std::abs(flow_controlled.network.measurements.pipeline_export_flow_kg_s -
                     pressure_boundary_step.measurements.pipeline_export_flow_kg_s) <
                1.0e-4,
            "pipeline-flow controller missed its calculated export-flow target");
    require(std::abs(flow_controlled.network.measurements.discharge_pressure_pa -
                     delivery_pressure_pa) < 2.0,
            "flow-controlled step did not preserve the delivery-pressure boundary");
    SteadySetpointOptions control_options;
    control_options.minimum_speed_rpm = 11990.0;
    control_options.maximum_speed_rpm = 12010.0;
    control_options.absolute_tolerance = 1.0e-3;
    control_options.relative_tolerance = 1.0e-10;
    control_options.maximum_iterations = 60;
    auto upstream_state = upstream_initial;
    auto downstream_state = downstream_initial;
    ControlledTwoPipeCompressorStep controlled;
    controlled = advance_controlled_two_pipe_compressor(
        SteadySetpointMode::DischargePressure,
        direct.measurements.discharge_pressure_pa,
        control_options,
        pipe, upstream_state, 500000.0,
        pipe, downstream_state, 5.0,
        compressor, coupling_options);
    require(controlled.control.status == SteadySetpointStatus::Converged,
            "speed controller did not converge against the downstream pipe tap");
    require(controlled.network.valid, controlled.network.message.c_str());
    require(std::abs(controlled.control.speed_rpm - 12000.0) < 1.0,
            "pressure controller did not recover the target compressor speed");
    require(std::abs(pipe_mass_balance_defect(
                pipe, upstream_state, controlled.network.upstream)) < 1.0e-6,
            "upstream pipe linepack balance failed on controlled step");
    require(std::abs(pipe_mass_balance_defect(
                pipe, downstream_state, controlled.network.downstream)) < 1.0e-6,
            "downstream pipe linepack balance failed on controlled step");


    const auto suction_controlled = advance_controlled_two_pipe_compressor(
        SteadySetpointMode::SuctionPressure,
        direct.measurements.suction_pressure_pa,
        control_options,
        pipe, upstream_initial, 500000.0,
        pipe, downstream_initial, 5.0,
        compressor, coupling_options);
    require(suction_controlled.control.status == SteadySetpointStatus::Converged,
            "suction-pressure controller did not converge against the upstream pipe tap");
    require(suction_controlled.network.valid,
            suction_controlled.network.message.c_str());
    require(std::abs(suction_controlled.control.speed_rpm - 12000.0) < 1.0,
            "suction-pressure controller did not recover the target compressor speed");
    require(std::abs(suction_controlled.network.measurements.suction_pressure_pa -
                     direct.measurements.suction_pressure_pa) < 1.0e-3,
            "suction-pressure controller missed its pipe-tap setpoint");
    auto final_step = controlled.network;
    upstream_state = final_step.upstream;
    downstream_state = final_step.downstream;
    for (int timestep = 1; timestep < 4; ++timestep) {
        final_step = advance_two_pipe_compressor(
            pipe, upstream_state, 500000.0,
            pipe, downstream_state, 5.5,
            compressor, 12000.0, coupling_options);
        require(final_step.valid, final_step.message.c_str());
        require(std::abs(final_step.compressor_node_mass_residual_kg_s) < 1.0e-6,
                "compressor-node mass balance failed after export-flow change");
        require(std::abs(pipe_mass_balance_defect(
                    pipe, upstream_state, final_step.upstream)) < 1.0e-6,
                "upstream pipe linepack balance failed after export-flow change");
        require(std::abs(pipe_mass_balance_defect(
                    pipe, downstream_state, final_step.downstream)) < 1.0e-6,
                "downstream pipe linepack balance failed after export-flow change");
        upstream_state = final_step.upstream;
        downstream_state = final_step.downstream;
    }

    const auto& result = final_step;
    require(result.compressor.valid, "compressor component did not converge");
    require(result.upstream.pressure_pa.size() == 3,
            "upstream pressure state has wrong size");
    require(result.downstream.pressure_pa.size() == 3,
            "downstream pressure state has wrong size");
    require(result.measurements.suction_pressure_pa > 0.0,
            "invalid compressor suction tap pressure");
    require(result.measurements.discharge_pressure_pa > 0.0,
            "invalid downstream delivery tap pressure");
    require(result.measurements.compressor_shaft_power_w > 0.0,
            "invalid compressor shaft power");
    require(result.measurements.pipeline_export_flow_kg_s == 5.5,
            "pipe export flow was not preserved through both legs");
    require(std::abs(result.compressor_node_mass_residual_kg_s) < 1.0e-6,
            "compressor-node mass balance did not close");
    require(std::abs(
                result.compressor.compressor_throughflow_kg_s -
                result.downstream.inlet_flow_kg_s) < 1.0e-6,
            "compressor flow does not match downstream pipe inlet flow");
    require(std::abs(
                result.compressor.compressor_throughflow_kg_s -
                result.measurements.pipeline_export_flow_kg_s) > 1.0e-6,
            "transient test did not distinguish compressor flow from export flow");
    require(std::abs(
                result.downstream.pressure_pa.front() -
                result.compressor.outlet_static_state->p) < 1.0e-8,
            "compressor discharge pressure was not applied to downstream pipe");
}

} // namespace

int main()
{
    test_two_pipe_compressor_step();
    return 0;
}
