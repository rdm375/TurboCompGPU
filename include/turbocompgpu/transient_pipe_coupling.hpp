#pragma once

#include <turbocompgpu/pipeline_compressor.hpp>

#include <pipe_sim/transient_step.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace turbocompgpu {

// One pressure/flow state for the isothermal pipe_sim discretization.
struct PipeSegmentState {
    std::vector<double> pressure_pa; // n + 1 pressure nodes
    std::vector<double> flow_kg_s;   // n staggered flow nodes
    double inlet_flow_kg_s = 0.0;
    double outlet_flow_kg_s = 0.0;   // prior imposed boundary flow
};

struct PipeSegmentModel {
    int cells = 0;
    pipe_sim::TransientParameters parameters{};
    pipe_sim::NewtonOptions solver_options{};
};

struct PipeSegmentStep {
    bool valid = false;
    std::string message;
    PipeSegmentState state;
    int newton_iterations = 0;
};

// Advances one isothermal pipe segment by one implicit theta-method step.
// The pipe solver owns its equations; this adapter only translates the
// segment state and boundary data into pipe_sim's public timestep API.
inline PipeSegmentStep advance_pipe_segment(
    const PipeSegmentModel& model,
    const PipeSegmentState& old_state,
    double inlet_pressure_pa,
    double outlet_flow_kg_s)
{
    const int n = model.cells;
    if (n < 2 || old_state.pressure_pa.size() !=
                     static_cast<std::size_t>(n + 1) ||
        old_state.flow_kg_s.size() != static_cast<std::size_t>(n) ||
        !std::isfinite(old_state.inlet_flow_kg_s) ||
        !std::isfinite(old_state.outlet_flow_kg_s) ||
        !std::isfinite(inlet_pressure_pa) || inlet_pressure_pa <= 0.0 ||
        !std::isfinite(outlet_flow_kg_s)) {
        throw std::invalid_argument("Invalid pipe segment step input");
    }

    std::vector<double> state(static_cast<std::size_t>(2 * n + 2));
    pipe_sim::NewtonWorkspace workspace(n);
    const pipe_sim::TransientBoundary boundary{
        old_state.inlet_flow_kg_s,
        old_state.outlet_flow_kg_s,
        outlet_flow_kg_s,
        inlet_pressure_pa
    };
    const auto step = pipe_sim::transient_step(
        n, old_state.pressure_pa, old_state.flow_kg_s, boundary,
        model.parameters, model.solver_options, state, workspace);

    PipeSegmentStep result;
    result.newton_iterations = step.iterations;
    if (step.info != 0) {
        result.message = "Pipe segment Newton step failed with info=" +
                         std::to_string(step.info);
        return result;
    }

    result.state.pressure_pa.resize(static_cast<std::size_t>(n + 1));
    result.state.flow_kg_s.resize(static_cast<std::size_t>(n));
    result.state.pressure_pa[0] = state[0];
    result.state.inlet_flow_kg_s = state[1];
    result.state.outlet_flow_kg_s = outlet_flow_kg_s;
    for (int i = 0; i < n; ++i) {
        result.state.flow_kg_s[static_cast<std::size_t>(i)] =
            state[static_cast<std::size_t>(2 * i + 2)];
        result.state.pressure_pa[static_cast<std::size_t>(i + 1)] =
            state[static_cast<std::size_t>(2 * i + 3)];
    }
    result.valid = true;
    result.message = "Pipe segment step converged";
    return result;
}

struct TwoPipeCompressorStep {
    bool valid = false;
    std::string message;
    PipeSegmentState upstream;
    PipeSegmentState downstream;
    CompressorComponentResult compressor;
    SteadyPlantMeasurements measurements;
    double compressor_node_mass_residual_kg_s = 0.0;
};

struct TwoPipeCompressorStepOptions {
    double minimum_compressor_throughflow_kg_s = 0.0;
    double maximum_compressor_throughflow_kg_s = 0.0;
    double absolute_mass_tolerance_kg_s = 1.0e-8;
    double relative_mass_tolerance = 1.0e-7;
    int maximum_flow_iterations = 60;
};

struct DeliveryPressureBoundaryOptions {
    TwoPipeCompressorStepOptions compressor_coupling;
    double minimum_export_flow_kg_s = 0.0;
    double maximum_export_flow_kg_s = 0.0;
    double absolute_pressure_tolerance_pa = 1.0;
    double relative_pressure_tolerance = 1.0e-8;
    int maximum_export_flow_iterations = 60;
};

struct TwoPipeRecycleCouplingOptions {
    double minimum_fresh_flow_kg_s = 0.0;
    double maximum_fresh_flow_kg_s = 0.0;
    double minimum_compressor_throughflow_kg_s = 0.0;
    double maximum_compressor_throughflow_kg_s = 0.0;
    double absolute_mass_tolerance_kg_s = 1.0e-8;
    double relative_mass_tolerance = 1.0e-7;
    int maximum_flow_iterations = 60;
    RecycleSteadyOptions recycle_loop;
};

struct TwoPipeRecycleCompressorStep {
    bool valid = false;
    std::string message;
    PipeSegmentState upstream;
    PipeSegmentState downstream;
    CompressorComponentResult compressor;
    SteadyPlantMeasurements measurements;
    int recycle_iterations = 0;
    double fresh_flow_kg_s = 0.0;
    double compressor_throughflow_kg_s = 0.0;
    double recycle_flow_kg_s = 0.0;
    double recycle_mass_residual_kg_s = 0.0;
    double recycle_cooler_heat_removed_w = 0.0;
    double suction_mixer_enthalpy_residual_w = 0.0;
    double fresh_flow_closure_residual_kg_s = 0.0;
    double compressor_node_mass_residual_kg_s = 0.0;
};

struct RecycleDeliveryPressureBoundaryOptions {
    TwoPipeRecycleCouplingOptions recycle_coupling;
    double minimum_export_flow_kg_s = 0.0;
    double maximum_export_flow_kg_s = 0.0;
    double absolute_pressure_tolerance_pa = 1.0;
    double relative_pressure_tolerance = 1.0e-8;
    int maximum_export_flow_iterations = 60;
};

struct ControlledTwoPipeRecycleCompressorStep {
    SteadySetpointResult control;
    TwoPipeRecycleCompressorStep network;
};

// Advances two pipe segments connected by a compressor for one time step.
// The upstream pipe has a fixed source pressure and trial compressor
// through-flow. The downstream pipe has the compressor discharge pressure
// and prescribed net export flow. A scalar inner solve enforces equality
// between compressor through-flow and the downstream pipe's computed inlet
// flow, so changing linepack is accounted for even with recycle closed.
//
// pipe_sim is isothermal. Each segment therefore uses its configured fixed
// temperature and compressibility factor; this adapter does not transport
// compressor discharge enthalpy into a non-isothermal pipe model.
inline TwoPipeCompressorStep advance_two_pipe_compressor(
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double pipeline_export_flow_kg_s,
    const CompressorComponent& compressor,
    double speed_rpm,
    const TwoPipeCompressorStepOptions& options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    if (!(std::isfinite(pipeline_export_flow_kg_s) &&
          pipeline_export_flow_kg_s > 0.0 &&
          std::isfinite(speed_rpm) && speed_rpm > 0.0 &&
          std::isfinite(upstream_model.parameters.temperature) &&
          upstream_model.parameters.temperature > 0.0 &&
          std::isfinite(options.minimum_compressor_throughflow_kg_s) &&
          options.minimum_compressor_throughflow_kg_s > 0.0 &&
          std::isfinite(options.maximum_compressor_throughflow_kg_s) &&
          options.maximum_compressor_throughflow_kg_s >
              options.minimum_compressor_throughflow_kg_s &&
          std::isfinite(options.absolute_mass_tolerance_kg_s) &&
          options.absolute_mass_tolerance_kg_s > 0.0 &&
          std::isfinite(options.relative_mass_tolerance) &&
          options.relative_mass_tolerance > 0.0 &&
          options.maximum_flow_iterations > 0)) {
        throw std::invalid_argument("Invalid two-pipe compressor request");
    }

    struct Trial {
        PipeSegmentStep upstream;
        PipeSegmentStep downstream;
        CompressorComponentResult compressor;
        double flow_residual_kg_s = 0.0;
        double suction_pressure_pa = 0.0;
    };
    auto evaluate_flow = [&](double throughflow) {
        Trial trial;
        trial.upstream = advance_pipe_segment(
            upstream_model, upstream_old, source_pressure_pa, throughflow);
        if (!trial.upstream.valid) {
            throw std::runtime_error("Upstream pipe: " + trial.upstream.message);
        }

        trial.suction_pressure_pa = trial.upstream.state.pressure_pa.back();
        trial.compressor = compressor.evaluate({
            {upstream_model.parameters.temperature, trial.suction_pressure_pa},
            throughflow,
            speed_rpm,
            guide_vane_angle_deg
        });
        if (!trial.compressor.valid ||
            !trial.compressor.outlet_static_state) {
            throw std::runtime_error("Compressor: " + trial.compressor.message);
        }

        trial.downstream = advance_pipe_segment(
            downstream_model, downstream_old,
            trial.compressor.outlet_static_state->p,
            pipeline_export_flow_kg_s);
        if (!trial.downstream.valid) {
            throw std::runtime_error("Downstream pipe: " + trial.downstream.message);
        }
        trial.flow_residual_kg_s =
            trial.downstream.state.inlet_flow_kg_s - throughflow;
        return trial;
    };

    double low = options.minimum_compressor_throughflow_kg_s;
    double high = options.maximum_compressor_throughflow_kg_s;
    auto flow_tolerance = [&](double flow) {
        return options.absolute_mass_tolerance_kg_s +
            options.relative_mass_tolerance * std::max(1.0, std::abs(flow));
    };

    const double prior_flow_guess = std::clamp(
        downstream_old.inlet_flow_kg_s, low, high);
    auto guess_trial = evaluate_flow(prior_flow_guess);
    std::optional<Trial> selected;
    if (std::abs(guess_trial.flow_residual_kg_s) <=
        flow_tolerance(prior_flow_guess)) {
        selected = std::move(guess_trial);
    }

    double bracket_low = 0.0;
    double bracket_high = 0.0;
    std::optional<Trial> bracket_low_trial;
    std::optional<Trial> bracket_high_trial;
    if (!selected) {
        constexpr int scan_intervals = 32;
        std::vector<std::pair<double, Trial>> samples;
        samples.reserve(scan_intervals + 1);
        for (int i = 0; i <= scan_intervals; ++i) {
            const double flow = low + (high - low) *
                static_cast<double>(i) / static_cast<double>(scan_intervals);
            auto trial = flow == prior_flow_guess
                ? std::move(guess_trial)
                : evaluate_flow(flow);
            if (std::abs(trial.flow_residual_kg_s) <= flow_tolerance(flow)) {
                selected = std::move(trial);
                break;
            }
            samples.emplace_back(flow, std::move(trial));
        }

        double best_bracket_distance =
            std::numeric_limits<double>::infinity();
        if (!selected) {
            for (std::size_t i = 1; i < samples.size(); ++i) {
                const auto& left = samples[i - 1];
                const auto& right = samples[i];
                if (std::signbit(left.second.flow_residual_kg_s) ==
                    std::signbit(right.second.flow_residual_kg_s)) {
                    continue;
                }
                const double distance = std::abs(
                    0.5 * (left.first + right.first) - prior_flow_guess);
                if (distance < best_bracket_distance) {
                    best_bracket_distance = distance;
                    bracket_low = left.first;
                    bracket_high = right.first;
                    bracket_low_trial = left.second;
                    bracket_high_trial = right.second;
                }
            }
        }
    }

    if (!selected && (!bracket_low_trial || !bracket_high_trial)) {
        TwoPipeCompressorStep failed;
        failed.message =
            "Compressor-flow closure has no root in the configured flow bounds";
        return failed;
    }

    if (!selected) {
        double local_low = bracket_low;
        double local_high = bracket_high;
        auto local_low_trial = std::move(*bracket_low_trial);
        auto local_high_trial = std::move(*bracket_high_trial);
        bool converged = false;
        for (int iteration = 0;
             iteration < options.maximum_flow_iterations;
             ++iteration) {
            const double mid = local_low + 0.5 * (local_high - local_low);
            auto mid_trial = evaluate_flow(mid);
            if (std::abs(mid_trial.flow_residual_kg_s) <= flow_tolerance(mid)) {
                selected = std::move(mid_trial);
                converged = true;
                break;
            }
            if (std::signbit(mid_trial.flow_residual_kg_s) ==
                std::signbit(local_low_trial.flow_residual_kg_s)) {
                local_low = mid;
                local_low_trial = std::move(mid_trial);
            } else {
                local_high = mid;
                local_high_trial = std::move(mid_trial);
            }
        }
        if (!converged) {
            TwoPipeCompressorStep failed;
            failed.message = "Compressor-flow closure iteration limit reached";
            return failed;
        }
    }

    TwoPipeCompressorStep result;
    result.valid = true;
    result.message = "Two-pipe compressor step converged";
    result.compressor_node_mass_residual_kg_s =
        selected->flow_residual_kg_s;
    result.measurements = {
        selected->suction_pressure_pa,
        selected->downstream.state.pressure_pa.back(),
        pipeline_export_flow_kg_s,
        selected->compressor.shaft_power_w
    };
    result.upstream = std::move(selected->upstream.state);
    result.downstream = std::move(selected->downstream.state);
    result.compressor = std::move(selected->compressor);
    return result;
}

// Advances a transient pipe/compressor network with a discharge-to-suction
// recycle branch. The upstream pipe carries fresh feed, while the compressor
// sees fresh feed plus recycle. The downstream pipe's computed inlet flow is
// closed against fresh feed, preserving both compressor-node mass balance and
// pipe linepack accounting. The recycle valve and return-state model must be
// supplied from machine/system calibration.
template<class RecycleReturnModel>
inline TwoPipeRecycleCompressorStep advance_two_pipe_compressor_with_recycle(
    const thermomodel::Model& eos,
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double pipeline_export_flow_kg_s,
    const CompressorComponent& compressor,
    double speed_rpm,
    double recycle_valve_opening,
    const RecycleValveMap& recycle_valve,
    RecycleReturnModel&& recycle_return_model,
    const TwoPipeRecycleCouplingOptions& options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    if (!(std::isfinite(pipeline_export_flow_kg_s) &&
          pipeline_export_flow_kg_s > 0.0 &&
          std::isfinite(speed_rpm) && speed_rpm > 0.0 &&
          std::isfinite(recycle_valve_opening) &&
          recycle_valve_opening >= 0.0 && recycle_valve_opening <= 1.0 &&
          std::isfinite(options.minimum_fresh_flow_kg_s) &&
          options.minimum_fresh_flow_kg_s > 0.0 &&
          std::isfinite(options.maximum_fresh_flow_kg_s) &&
          options.maximum_fresh_flow_kg_s > options.minimum_fresh_flow_kg_s &&
          std::isfinite(options.minimum_compressor_throughflow_kg_s) &&
          options.minimum_compressor_throughflow_kg_s > 0.0 &&
          std::isfinite(options.maximum_compressor_throughflow_kg_s) &&
          options.maximum_compressor_throughflow_kg_s >
              options.minimum_compressor_throughflow_kg_s &&
          std::isfinite(options.absolute_mass_tolerance_kg_s) &&
          options.absolute_mass_tolerance_kg_s > 0.0 &&
          std::isfinite(options.relative_mass_tolerance) &&
          options.relative_mass_tolerance > 0.0 &&
          options.maximum_flow_iterations > 0)) {
        throw std::invalid_argument("Invalid two-pipe recycle request");
    }
    recycle_valve.validate();

    struct Trial {
        PipeSegmentStep upstream;
        PipeSegmentStep downstream;
        RecycleSteadyResult recycle;
        double fresh_flow_kg_s = 0.0;
        double closure_residual_kg_s = 0.0;
        double suction_pressure_pa = 0.0;
    };
    auto evaluate_fresh_flow = [&](double fresh_flow) {
        Trial trial;
        trial.fresh_flow_kg_s = fresh_flow;
        trial.upstream = advance_pipe_segment(
            upstream_model, upstream_old, source_pressure_pa, fresh_flow);
        if (!trial.upstream.valid) {
            throw std::runtime_error("Upstream pipe: " + trial.upstream.message);
        }
        trial.suction_pressure_pa = trial.upstream.state.pressure_pa.back();
        const CompressorPipePort fresh_inlet{
            upstream_model.parameters.temperature,
            trial.suction_pressure_pa
        };

        auto evaluate_compressor = [&](const thermomodel::State& mixed,
                                       double throughflow) {
            if (throughflow < options.minimum_compressor_throughflow_kg_s ||
                throughflow > options.maximum_compressor_throughflow_kg_s) {
                CompressorComponentResult invalid;
                invalid.message =
                    "Recycle loop compressor through-flow is outside configured bounds";
                return invalid;
            }
            return compressor.evaluate({
                {mixed.T, mixed.p}, throughflow, speed_rpm,
                guide_vane_angle_deg
            });
        };
        trial.recycle = solve_steady_recycle(
            eos, fresh_inlet, fresh_flow, recycle_valve_opening,
            recycle_valve, evaluate_compressor,
            recycle_return_model, options.recycle_loop);
        if (!trial.recycle.converged ||
            !trial.recycle.compressor.outlet_static_state) {
            throw std::runtime_error("Recycle loop: " + trial.recycle.message);
        }
        trial.downstream = advance_pipe_segment(
            downstream_model, downstream_old,
            trial.recycle.compressor.outlet_static_state->p,
            pipeline_export_flow_kg_s);
        if (!trial.downstream.valid) {
            throw std::runtime_error("Downstream pipe: " +
                                     trial.downstream.message);
        }
        trial.closure_residual_kg_s =
            trial.downstream.state.inlet_flow_kg_s - fresh_flow;
        return trial;
    };

    const double low = options.minimum_fresh_flow_kg_s;
    const double high = options.maximum_fresh_flow_kg_s;
    auto tolerance = [&](double flow) {
        return options.absolute_mass_tolerance_kg_s +
            options.relative_mass_tolerance * std::max(1.0, std::abs(flow));
    };
    const double prior_guess = std::clamp(
        downstream_old.inlet_flow_kg_s, low, high);
    auto guess = evaluate_fresh_flow(prior_guess);
    std::optional<Trial> selected;
    if (std::abs(guess.closure_residual_kg_s) <= tolerance(prior_guess)) {
        selected = std::move(guess);
    }

    std::vector<std::pair<double, Trial>> samples;
    if (!selected) {
        constexpr int scan_intervals = 32;
        samples.reserve(scan_intervals + 1);
        for (int i = 0; i <= scan_intervals; ++i) {
            const double flow = low + (high - low) *
                static_cast<double>(i) / static_cast<double>(scan_intervals);
            auto trial = flow == prior_guess
                ? std::move(guess) : evaluate_fresh_flow(flow);
            if (std::abs(trial.closure_residual_kg_s) <= tolerance(flow)) {
                selected = std::move(trial);
                break;
            }
            samples.emplace_back(flow, std::move(trial));
        }
    }

    double bracket_low = 0.0;
    double bracket_high = 0.0;
    std::optional<Trial> bracket_low_trial;
    std::optional<Trial> bracket_high_trial;
    double best_distance = std::numeric_limits<double>::infinity();
    if (!selected) {
        for (std::size_t i = 1; i < samples.size(); ++i) {
            const auto& left = samples[i - 1];
            const auto& right = samples[i];
            if (std::signbit(left.second.closure_residual_kg_s) ==
                std::signbit(right.second.closure_residual_kg_s)) {
                continue;
            }
            const double distance = std::abs(
                0.5 * (left.first + right.first) - prior_guess);
            if (distance < best_distance) {
                best_distance = distance;
                bracket_low = left.first;
                bracket_high = right.first;
                bracket_low_trial = left.second;
                bracket_high_trial = right.second;
            }
        }
    }
    if (!selected && (!bracket_low_trial || !bracket_high_trial)) {
        TwoPipeRecycleCompressorStep failed;
        failed.message =
            "Fresh-flow closure has no root in the configured flow bounds";
        return failed;
    }
    if (!selected) {
        double local_low = bracket_low;
        double local_high = bracket_high;
        auto local_low_trial = std::move(*bracket_low_trial);
        bool converged = false;
        for (int i = 0; i < options.maximum_flow_iterations; ++i) {
            const double mid = local_low + 0.5 * (local_high - local_low);
            auto mid_trial = evaluate_fresh_flow(mid);
            if (std::abs(mid_trial.closure_residual_kg_s) <= tolerance(mid)) {
                selected = std::move(mid_trial);
                converged = true;
                break;
            }
            if (std::signbit(mid_trial.closure_residual_kg_s) ==
                std::signbit(local_low_trial.closure_residual_kg_s)) {
                local_low = mid;
                local_low_trial = std::move(mid_trial);
            } else {
                local_high = mid;
            }
        }
        if (!converged) {
            TwoPipeRecycleCompressorStep failed;
            failed.message = "Fresh-flow closure iteration limit reached";
            return failed;
        }
    }

    TwoPipeRecycleCompressorStep result;
    result.valid = true;
    result.message = "Two-pipe recycle compressor step converged";
    result.fresh_flow_closure_residual_kg_s =
        selected->closure_residual_kg_s;
    result.compressor_node_mass_residual_kg_s =
        selected->recycle.compressor_throughflow_kg_s -
        selected->recycle.recycle_flow_kg_s -
        selected->downstream.state.inlet_flow_kg_s;
    result.measurements = {
        selected->suction_pressure_pa,
        selected->downstream.state.pressure_pa.back(),
        pipeline_export_flow_kg_s,
        selected->recycle.compressor.shaft_power_w
    };
    result.recycle_iterations = selected->recycle.iterations;
    result.fresh_flow_kg_s = selected->fresh_flow_kg_s;
    result.compressor_throughflow_kg_s =
        selected->recycle.compressor_throughflow_kg_s;
    result.recycle_flow_kg_s = selected->recycle.recycle_flow_kg_s;
    result.recycle_mass_residual_kg_s =
        selected->recycle.recycle_mass_residual_kg_s;
    result.recycle_cooler_heat_removed_w =
        selected->recycle.recycle_cooler_heat_removed_w;
    result.suction_mixer_enthalpy_residual_w =
        selected->recycle.suction_mixer_enthalpy_residual_w;
    result.upstream = std::move(selected->upstream.state);
    result.downstream = std::move(selected->downstream.state);
    result.compressor = std::move(selected->recycle.compressor);
    return result;
}

// Solves the downstream boundary flow required to hold a delivery pressure.
// This complements the prescribed-export-flow form above: export is now a
// network result and can be used meaningfully by a pipeline-flow controller.
inline TwoPipeCompressorStep advance_two_pipe_compressor_at_delivery_pressure(
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double delivery_pressure_pa,
    const CompressorComponent& compressor,
    double speed_rpm,
    const DeliveryPressureBoundaryOptions& options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    if (!(std::isfinite(delivery_pressure_pa) && delivery_pressure_pa > 0.0 &&
          std::isfinite(options.minimum_export_flow_kg_s) &&
          options.minimum_export_flow_kg_s > 0.0 &&
          std::isfinite(options.maximum_export_flow_kg_s) &&
          options.maximum_export_flow_kg_s >
              options.minimum_export_flow_kg_s &&
          std::isfinite(options.absolute_pressure_tolerance_pa) &&
          options.absolute_pressure_tolerance_pa > 0.0 &&
          std::isfinite(options.relative_pressure_tolerance) &&
          options.relative_pressure_tolerance > 0.0 &&
          options.maximum_export_flow_iterations > 0)) {
        throw std::invalid_argument("Invalid delivery-pressure boundary request");
    }

    auto pressure_tolerance = [&](double pressure) {
        return options.absolute_pressure_tolerance_pa +
            options.relative_pressure_tolerance * std::abs(pressure);
    };
    auto evaluate_export = [&](double export_flow) {
        return advance_two_pipe_compressor(
            upstream_model, upstream_old, source_pressure_pa,
            downstream_model, downstream_old, export_flow,
            compressor, speed_rpm, options.compressor_coupling,
            guide_vane_angle_deg);
    };

    struct Sample {
        double export_flow = 0.0;
        double pressure_residual_pa = 0.0;
        TwoPipeCompressorStep step;
    };
    constexpr int scan_intervals = 32;
    std::vector<Sample> samples;
    samples.reserve(scan_intervals + 1);
    std::optional<TwoPipeCompressorStep> selected;
    const double prior_export = std::clamp(
        downstream_old.outlet_flow_kg_s,
        options.minimum_export_flow_kg_s,
        options.maximum_export_flow_kg_s);
    for (int i = 0; i <= scan_intervals; ++i) {
        const double export_flow = options.minimum_export_flow_kg_s +
            (options.maximum_export_flow_kg_s -
             options.minimum_export_flow_kg_s) *
            static_cast<double>(i) / static_cast<double>(scan_intervals);
        auto step = evaluate_export(export_flow);
        if (!step.valid) continue;
        const double residual =
            step.measurements.discharge_pressure_pa - delivery_pressure_pa;
        if (std::abs(residual) <= pressure_tolerance(delivery_pressure_pa)) {
            if (!selected ||
                std::abs(export_flow - prior_export) <
                    std::abs(selected->measurements.pipeline_export_flow_kg_s -
                             prior_export)) {
                selected = std::move(step);
            }
            continue;
        }
        samples.push_back({export_flow, residual, std::move(step)});
    }

    double bracket_low = 0.0;
    double bracket_high = 0.0;
    std::optional<Sample> low_sample;
    std::optional<Sample> high_sample;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < samples.size(); ++i) {
        if (std::signbit(samples[i - 1].pressure_residual_pa) ==
            std::signbit(samples[i].pressure_residual_pa)) {
            continue;
        }
        const double distance = std::abs(
            0.5 * (samples[i - 1].export_flow + samples[i].export_flow) -
            prior_export);
        if (distance < best_distance) {
            best_distance = distance;
            bracket_low = samples[i - 1].export_flow;
            bracket_high = samples[i].export_flow;
            low_sample = samples[i - 1];
            high_sample = samples[i];
        }
    }

    if (!selected && (!low_sample || !high_sample)) {
        TwoPipeCompressorStep failed;
        failed.message =
            "Delivery pressure is not bracketed by the configured export-flow bounds";
        return failed;
    }

    if (!selected) {
        double low = bracket_low;
        double high = bracket_high;
        double low_residual = low_sample->pressure_residual_pa;
        for (int iteration = 0;
             iteration < options.maximum_export_flow_iterations;
             ++iteration) {
            const double mid = low + 0.5 * (high - low);
            auto step = evaluate_export(mid);
            if (!step.valid) {
                TwoPipeCompressorStep failed;
                failed.message = "Pipe-compressor trial failed during delivery-pressure solve: " +
                                step.message;
                return failed;
            }
            const double residual =
                step.measurements.discharge_pressure_pa - delivery_pressure_pa;
            if (std::abs(residual) <= pressure_tolerance(delivery_pressure_pa)) {
                selected = std::move(step);
                break;
            }
            if (std::signbit(residual) == std::signbit(low_residual)) {
                low = mid;
                low_residual = residual;
            } else {
                high = mid;
            }
        }
        if (!selected) {
            TwoPipeCompressorStep failed;
            failed.message = "Delivery-pressure boundary iteration limit reached";
            return failed;
        }
    }
    return std::move(*selected);
}

// Recycle-aware form of the delivery-pressure boundary. Export flow is
// calculated from the requested downstream pressure; the recycle-aware
// inner solve still closes compressor through-flow as fresh feed plus recycle.
template<class RecycleReturnModel>
inline TwoPipeRecycleCompressorStep
advance_two_pipe_compressor_with_recycle_at_delivery_pressure(
    const thermomodel::Model& eos,
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double delivery_pressure_pa,
    const CompressorComponent& compressor,
    double speed_rpm,
    double recycle_valve_opening,
    const RecycleValveMap& recycle_valve,
    RecycleReturnModel&& recycle_return_model,
    const RecycleDeliveryPressureBoundaryOptions& options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    if (!(std::isfinite(delivery_pressure_pa) && delivery_pressure_pa > 0.0 &&
          std::isfinite(options.minimum_export_flow_kg_s) &&
          options.minimum_export_flow_kg_s > 0.0 &&
          std::isfinite(options.maximum_export_flow_kg_s) &&
          options.maximum_export_flow_kg_s > options.minimum_export_flow_kg_s &&
          std::isfinite(options.absolute_pressure_tolerance_pa) &&
          options.absolute_pressure_tolerance_pa > 0.0 &&
          std::isfinite(options.relative_pressure_tolerance) &&
          options.relative_pressure_tolerance > 0.0 &&
          options.maximum_export_flow_iterations > 0)) {
        throw std::invalid_argument(
            "Invalid recycle delivery-pressure boundary request");
    }

    auto tolerance = [&](double pressure) {
        return options.absolute_pressure_tolerance_pa +
            options.relative_pressure_tolerance * std::abs(pressure);
    };
    auto evaluate_export = [&](double export_flow) {
        return advance_two_pipe_compressor_with_recycle(
            eos, upstream_model, upstream_old, source_pressure_pa,
            downstream_model, downstream_old, export_flow,
            compressor, speed_rpm, recycle_valve_opening,
            recycle_valve, recycle_return_model,
            options.recycle_coupling, guide_vane_angle_deg);
    };

    const double low = options.minimum_export_flow_kg_s;
    const double high = options.maximum_export_flow_kg_s;
    const double prior_export = std::clamp(
        downstream_old.outlet_flow_kg_s, low, high);
    auto prior = evaluate_export(prior_export);
    if (prior.valid && std::abs(
            prior.measurements.discharge_pressure_pa - delivery_pressure_pa) <=
            tolerance(delivery_pressure_pa)) {
        return prior;
    }

    struct Sample {
        double flow = 0.0;
        double residual_pa = 0.0;
        TwoPipeRecycleCompressorStep step;
    };
    std::vector<Sample> samples;
    samples.reserve(34);
    if (prior.valid) {
        samples.push_back({
            prior_export,
            prior.measurements.discharge_pressure_pa - delivery_pressure_pa,
            std::move(prior)
        });
    }
    constexpr int scan_intervals = 32;
    for (int i = 0; i <= scan_intervals; ++i) {
        const double flow = low + (high - low) *
            static_cast<double>(i) / static_cast<double>(scan_intervals);
        if (flow == prior_export) continue;
        auto step = evaluate_export(flow);
        if (!step.valid) continue;
        const double residual =
            step.measurements.discharge_pressure_pa - delivery_pressure_pa;
        if (std::abs(residual) <= tolerance(delivery_pressure_pa)) {
            return step;
        }
        samples.push_back({flow, residual, std::move(step)});
    }
    std::sort(samples.begin(), samples.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.flow < rhs.flow;
    });

    std::optional<Sample> bracket_low;
    std::optional<Sample> bracket_high;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < samples.size(); ++i) {
        if (std::signbit(samples[i - 1].residual_pa) ==
            std::signbit(samples[i].residual_pa)) {
            continue;
        }
        const double distance = std::abs(
            0.5 * (samples[i - 1].flow + samples[i].flow) - prior_export);
        if (distance < best_distance) {
            best_distance = distance;
            bracket_low = samples[i - 1];
            bracket_high = samples[i];
        }
    }
    if (!bracket_low || !bracket_high) {
        TwoPipeRecycleCompressorStep failed;
        failed.message =
            "Delivery pressure is not bracketed by recycle export-flow bounds";
        return failed;
    }

    double local_low = bracket_low->flow;
    double local_high = bracket_high->flow;
    double low_residual = bracket_low->residual_pa;
    for (int iteration = 0;
         iteration < options.maximum_export_flow_iterations;
         ++iteration) {
        const double mid = local_low + 0.5 * (local_high - local_low);
        auto step = evaluate_export(mid);
        if (!step.valid) {
            TwoPipeRecycleCompressorStep failed;
            failed.message = "Recycle network failed during delivery-pressure solve: " +
                             step.message;
            return failed;
        }
        const double residual =
            step.measurements.discharge_pressure_pa - delivery_pressure_pa;
        if (std::abs(residual) <= tolerance(delivery_pressure_pa)) {
            return step;
        }
        if (std::signbit(residual) == std::signbit(low_residual)) {
            local_low = mid;
            low_residual = residual;
        } else {
            local_high = mid;
        }
    }

    TwoPipeRecycleCompressorStep failed;
    failed.message = "Recycle delivery-pressure iteration limit reached";
    return failed;
}

struct ControlledTwoPipeCompressorStep {
    SteadySetpointResult control;
    TwoPipeCompressorStep network;
};

// Connects the A1 speed controller to the actual two-segment timestep. Each
// trial speed starts from the same old pipe states; only the final selected
// speed's pipe states are returned to the caller for committing.
inline ControlledTwoPipeCompressorStep advance_controlled_two_pipe_compressor(
    SteadySetpointMode mode,
    double setpoint,
    const SteadySetpointOptions& options,
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double pipeline_export_flow_kg_s,
    const CompressorComponent& compressor,
    const TwoPipeCompressorStepOptions& coupling_options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    auto evaluate = [&](double speed) {
        auto trial = advance_two_pipe_compressor(
            upstream_model, upstream_old, source_pressure_pa,
            downstream_model, downstream_old, pipeline_export_flow_kg_s,
            compressor, speed, coupling_options, guide_vane_angle_deg);
        if (!trial.valid) {
            throw std::runtime_error(
                "Coupled plant evaluation failed: " + trial.message);
        }
        return trial.measurements;
    };

    ControlledTwoPipeCompressorStep result;
    result.control = solve_steady_speed_setpoint(
        mode, setpoint, options, evaluate);
    result.network = advance_two_pipe_compressor(
        upstream_model, upstream_old, source_pressure_pa,
        downstream_model, downstream_old, pipeline_export_flow_kg_s,
        compressor, result.control.speed_rpm, coupling_options,
        guide_vane_angle_deg);
    return result;
}

// Pipeline-flow control against a fixed delivery-pressure boundary. Unlike
// the prescribed-export-flow variant, the controller observes the flow
// calculated by the network as its delivered-flow measurement.
inline ControlledTwoPipeCompressorStep
advance_flow_controlled_two_pipe_compressor_at_delivery_pressure(
    double pipeline_flow_setpoint_kg_s,
    const SteadySetpointOptions& control_options,
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double delivery_pressure_pa,
    const CompressorComponent& compressor,
    const DeliveryPressureBoundaryOptions& boundary_options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    auto evaluate = [&](double speed) {
        auto trial = advance_two_pipe_compressor_at_delivery_pressure(
            upstream_model, upstream_old, source_pressure_pa,
            downstream_model, downstream_old, delivery_pressure_pa,
            compressor, speed, boundary_options, guide_vane_angle_deg);
        if (!trial.valid) {
            throw std::runtime_error(
                "Coupled delivery-pressure evaluation failed: " + trial.message);
        }
        return trial.measurements;
    };

    ControlledTwoPipeCompressorStep result;
    result.control = solve_steady_speed_setpoint(
        SteadySetpointMode::PipelineExportFlow,
        pipeline_flow_setpoint_kg_s, control_options, evaluate);
    result.network = advance_two_pipe_compressor_at_delivery_pressure(
        upstream_model, upstream_old, source_pressure_pa,
        downstream_model, downstream_old, delivery_pressure_pa,
        compressor, result.control.speed_rpm, boundary_options,
        guide_vane_angle_deg);
    return result;
}
template<class RecycleReturnModel>
inline ControlledTwoPipeRecycleCompressorStep
advance_flow_controlled_two_pipe_compressor_with_recycle_at_delivery_pressure(
    double pipeline_flow_setpoint_kg_s,
    const SteadySetpointOptions& control_options,
    const thermomodel::Model& eos,
    const PipeSegmentModel& upstream_model,
    const PipeSegmentState& upstream_old,
    double source_pressure_pa,
    const PipeSegmentModel& downstream_model,
    const PipeSegmentState& downstream_old,
    double delivery_pressure_pa,
    const CompressorComponent& compressor,
    double recycle_valve_opening,
    const RecycleValveMap& recycle_valve,
    RecycleReturnModel&& recycle_return_model,
    const RecycleDeliveryPressureBoundaryOptions& boundary_options,
    std::optional<double> guide_vane_angle_deg = std::nullopt)
{
    auto evaluate = [&](double speed) {
        auto trial = advance_two_pipe_compressor_with_recycle_at_delivery_pressure(
            eos, upstream_model, upstream_old, source_pressure_pa,
            downstream_model, downstream_old, delivery_pressure_pa,
            compressor, speed, recycle_valve_opening, recycle_valve,
            recycle_return_model, boundary_options, guide_vane_angle_deg);
        if (!trial.valid) {
            throw std::runtime_error(
                "Recycle plant evaluation failed: " + trial.message);
        }
        return trial.measurements;
    };

    ControlledTwoPipeRecycleCompressorStep result;
    result.control = solve_steady_speed_setpoint(
        SteadySetpointMode::PipelineExportFlow,
        pipeline_flow_setpoint_kg_s, control_options, evaluate);
    result.network = advance_two_pipe_compressor_with_recycle_at_delivery_pressure(
        eos, upstream_model, upstream_old, source_pressure_pa,
        downstream_model, downstream_old, delivery_pressure_pa,
        compressor, result.control.speed_rpm, recycle_valve_opening,
        recycle_valve, recycle_return_model, boundary_options,
        guide_vane_angle_deg);
    return result;
}

} // namespace turbocompgpu
