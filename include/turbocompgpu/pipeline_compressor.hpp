#pragma once

#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>
#include <turbocompgpu/thermodynamic_solver.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace turbocompgpu {

// A0: Pipe ports carry static thermodynamic state. Positive flow is
// from the upstream pipe, through the compressor, to the downstream pipe.
struct CompressorPipePort {
    double static_temperature_k = 0.0;
    double static_pressure_pa = 0.0;
};

struct GuideVaneEffects {
    double inlet_prewhirl_m_s = 0.0;
    double entropy_generation_j_kg_k = 0.0;
};

struct GuideVaneOperatingDomain {
    double minimum_angle_deg = 0.0;
    double maximum_angle_deg = 0.0;
    double minimum_speed_rpm = 0.0;
    double maximum_speed_rpm = 0.0;
    double minimum_flow_kg_s = 0.0;
    double maximum_flow_kg_s = 0.0;
};

// A2: Calibration is explicit and bounded. The callback should come from
// machine data or a separately documented calibration; no vane physics is
// inferred from an arbitrary angle-to-prewhirl assumption.
struct GuideVaneMap {
    GuideVaneOperatingDomain domain;
    std::function<GuideVaneEffects(
        double angle_deg,
        double speed_rpm,
        double compressor_flow_kg_s,
        const thermomodel::State& inlet_state)> calibrated_effects;

    void validate() const
    {
        const auto& d = domain;
        if (!(std::isfinite(d.minimum_angle_deg) &&
              std::isfinite(d.maximum_angle_deg) &&
              d.maximum_angle_deg > d.minimum_angle_deg &&
              std::isfinite(d.minimum_speed_rpm) &&
              std::isfinite(d.maximum_speed_rpm) &&
              d.minimum_speed_rpm > 0.0 &&
              d.maximum_speed_rpm > d.minimum_speed_rpm &&
              std::isfinite(d.minimum_flow_kg_s) &&
              std::isfinite(d.maximum_flow_kg_s) &&
              d.minimum_flow_kg_s >= 0.0 &&
              d.maximum_flow_kg_s > d.minimum_flow_kg_s)) {
            throw std::invalid_argument("Invalid guide-vane calibration domain");
        }
    }

    GuideVaneEffects evaluate(
        double angle_deg,
        double speed_rpm,
        double compressor_flow_kg_s,
        const thermomodel::State& inlet_state) const
    {
        validate();
        if (!calibrated_effects ||
            !std::isfinite(angle_deg) ||
            !std::isfinite(speed_rpm) ||
            !std::isfinite(compressor_flow_kg_s) ||
            angle_deg < domain.minimum_angle_deg ||
            angle_deg > domain.maximum_angle_deg ||
            speed_rpm < domain.minimum_speed_rpm ||
            speed_rpm > domain.maximum_speed_rpm ||
            compressor_flow_kg_s < domain.minimum_flow_kg_s ||
            compressor_flow_kg_s > domain.maximum_flow_kg_s) {
            throw std::out_of_range(
                "Guide-vane query is outside its calibrated domain");
        }
        const auto effects = calibrated_effects(
            angle_deg, speed_rpm, compressor_flow_kg_s, inlet_state);
        if (!(std::isfinite(effects.inlet_prewhirl_m_s) &&
              std::isfinite(effects.entropy_generation_j_kg_k) &&
              effects.entropy_generation_j_kg_k >= 0.0)) {
            throw std::runtime_error(
                "Guide-vane calibration returned invalid effects");
        }
        return effects;
    }
};

struct CompressorComponentConfig {
    double slip_factor = 1.0;
    GuideVaneEffects fixed_vane_effects{};
    DiffuserBranch diffuser_branch = DiffuserBranch::Subsonic;
    NonlinearResidualScales residual_scales{1.0e5, 1.0e3, 1.0e5};
    NonlinearSolverOptions solver_options{};
    double validation_tolerance = 1.0e-8;
    double initial_outlet_temperature_k = 350.0;
    double initial_outlet_pressure_pa = 700000.0;
};

struct CompressorComponentInput {
    CompressorPipePort inlet;
    double compressor_throughflow_kg_s = 0.0;
    double speed_rpm = 0.0;
    std::optional<double> guide_vane_angle_deg;
};

struct CompressorComponentResult {
    bool valid = false;
    std::string message;
    double compressor_throughflow_kg_s = 0.0;
    double speed_rpm = 0.0;
    std::optional<double> guide_vane_angle_deg;
    GuideVaneEffects vane_effects;
    std::optional<thermomodel::State> outlet_static_state;
    std::optional<thermomodel::State> outlet_total_state;
    double shaft_power_w = std::numeric_limits<double>::quiet_NaN();
    CompleteStageResult stage;
};

class CompressorComponent {
public:
    CompressorComponent(
        const thermomodel::Model& model,
        const thermomodel::DifferentiableModel& derivatives,
        ImpellerSystemGeometry impeller_geometry,
        VanelessDiffuserGeometry diffuser_geometry,
        CompressorComponentConfig config,
        const GuideVaneMap* guide_vane_map = nullptr)
        : model_(model),
          derivatives_(derivatives),
          impeller_geometry_(std::move(impeller_geometry)),
          diffuser_geometry_(std::move(diffuser_geometry)),
          config_(std::move(config)),
          guide_vane_map_(guide_vane_map)
    {
        if (!(std::isfinite(config_.slip_factor) &&
              config_.slip_factor > 0.0 &&
              config_.slip_factor <= 1.0 &&
              std::isfinite(config_.validation_tolerance) &&
              config_.validation_tolerance > 0.0 &&
              std::isfinite(config_.residual_scales.specific_work) &&
              config_.residual_scales.specific_work > 0.0 &&
              std::isfinite(config_.residual_scales.specific_entropy) &&
              config_.residual_scales.specific_entropy > 0.0 &&
              std::isfinite(config_.residual_scales.pressure) &&
              config_.residual_scales.pressure > 0.0 &&
              std::isfinite(config_.fixed_vane_effects.inlet_prewhirl_m_s) &&
              std::isfinite(config_.fixed_vane_effects.entropy_generation_j_kg_k) &&
              config_.fixed_vane_effects.entropy_generation_j_kg_k >= 0.0 &&
              std::isfinite(config_.initial_outlet_temperature_k) &&
              config_.initial_outlet_temperature_k > 0.0 &&
              std::isfinite(config_.initial_outlet_pressure_pa) &&
              config_.initial_outlet_pressure_pa > 0.0)) {
            throw std::invalid_argument(
                "Invalid compressor-component configuration");
        }
    }

    CompressorComponentResult evaluate(
        const CompressorComponentInput& input) const
    {
        if (!(std::isfinite(input.inlet.static_temperature_k) &&
              input.inlet.static_temperature_k > 0.0 &&
              std::isfinite(input.inlet.static_pressure_pa) &&
              input.inlet.static_pressure_pa > 0.0 &&
              std::isfinite(input.compressor_throughflow_kg_s) &&
              input.compressor_throughflow_kg_s > 0.0 &&
              std::isfinite(input.speed_rpm) &&
              input.speed_rpm > 0.0)) {
            throw std::invalid_argument(
                "Invalid compressor-component input");
        }

        const auto inlet_state = turbocompgpu::evaluate(
            model_,
            input.inlet.static_temperature_k,
            input.inlet.static_pressure_pa);

        GuideVaneEffects vane = config_.fixed_vane_effects;
        if (input.guide_vane_angle_deg.has_value()) {
            if (guide_vane_map_ == nullptr) {
                throw std::invalid_argument(
                    "Guide-vane angle supplied without a calibrated map");
            }
            vane = guide_vane_map_->evaluate(
                *input.guide_vane_angle_deg,
                input.speed_rpm,
                input.compressor_throughflow_kg_s,
                inlet_state);
        }

        if (!(std::isfinite(vane.inlet_prewhirl_m_s) &&
              std::isfinite(vane.entropy_generation_j_kg_k) &&
              vane.entropy_generation_j_kg_k >= 0.0)) {
            throw std::invalid_argument(
                "Invalid fixed-vane effects");
        }

        const NonlinearOperatingPoint operating{
            input.inlet.static_temperature_k,
            input.inlet.static_pressure_pa,
            input.compressor_throughflow_kg_s,
            input.speed_rpm,
            config_.slip_factor,
            vane.inlet_prewhirl_m_s,
            vane.entropy_generation_j_kg_k
        };
        const OperatingPointRequest request{
            operating,
            {
                config_.initial_outlet_temperature_k,
                std::log(config_.initial_outlet_pressure_pa /
                         config_.residual_scales.pressure)
            }
        };

        CompressorComponentResult result;
        result.compressor_throughflow_kg_s =
            input.compressor_throughflow_kg_s;
        result.speed_rpm = input.speed_rpm;
        result.guide_vane_angle_deg = input.guide_vane_angle_deg;
        result.vane_effects = vane;
        result.stage = evaluate_complete_stage(
            model_,
            derivatives_,
            impeller_geometry_,
            diffuser_geometry_,
            request,
            config_.residual_scales,
            config_.solver_options,
            config_.diffuser_branch,
            config_.validation_tolerance);

        if (!result.stage.complete() ||
            !result.stage.stage ||
            !result.stage.stage->station3) {
            result.message = result.stage.message;
            return result;
        }

        const auto& outlet = *result.stage.stage->station3;
        result.outlet_static_state = outlet.state;
        const FlowStation outlet_station{
            outlet.state,
            {outlet.meridional_velocity, outlet.tangential_velocity}
        };
        result.outlet_total_state =
            evaluate_stagnation_state(model_, outlet_station).stagnation_state;
        result.shaft_power_w =
            result.stage.stage->impeller_shaft_power;
        result.valid = true;
        result.message = result.stage.message;
        return result;
    }

private:
    const thermomodel::Model& model_;
    const thermomodel::DifferentiableModel& derivatives_;
    ImpellerSystemGeometry impeller_geometry_;
    VanelessDiffuserGeometry diffuser_geometry_;
    CompressorComponentConfig config_;
    const GuideVaneMap* guide_vane_map_;
};

// A1: a pipe solver supplies these measured steady-state values after solving
// its network with the compressor at the requested speed.
struct SteadyPlantMeasurements {
    double suction_pressure_pa = 0.0;
    double discharge_pressure_pa = 0.0;
    double pipeline_export_flow_kg_s = 0.0;
    double compressor_shaft_power_w =
        std::numeric_limits<double>::quiet_NaN();
};

enum class SteadySetpointMode {
    SuctionPressure,
    DischargePressure,
    PipelineExportFlow,
    ManualSpeed
};

enum class SteadySetpointStatus {
    Converged,
    SpeedLimitReached,
    PowerLimitReached,
    TargetNotBracketed,
    IterationLimit
};

struct SteadySetpointOptions {
    double minimum_speed_rpm = 0.0;
    double maximum_speed_rpm = 0.0;
    double maximum_shaft_power_w =
        std::numeric_limits<double>::infinity();
    double absolute_tolerance = 0.0;
    double relative_tolerance = 1.0e-6;
    int maximum_iterations = 80;
};

struct SteadySetpointResult {
    SteadySetpointStatus status = SteadySetpointStatus::TargetNotBracketed;
    double speed_rpm = 0.0;
    double measured_value = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    int plant_evaluations = 0;
    SteadyPlantMeasurements measurements;
};

inline void validate_plant_measurements(
    const SteadyPlantMeasurements& measurements)
{
    if (!(std::isfinite(measurements.suction_pressure_pa) &&
          measurements.suction_pressure_pa > 0.0 &&
          std::isfinite(measurements.discharge_pressure_pa) &&
          measurements.discharge_pressure_pa > 0.0 &&
          std::isfinite(measurements.pipeline_export_flow_kg_s) &&
          measurements.pipeline_export_flow_kg_s >= 0.0)) {
        throw std::runtime_error(
            "Pipe/compressor evaluator returned invalid measurements");
    }
}

// The evaluator must solve the attached pipe/compressor/recycle network at
// the supplied speed. Power-limited speed search assumes shaft power is
// nondecreasing over the configured speed interval.
inline SteadySetpointResult solve_steady_speed_setpoint(
    SteadySetpointMode mode,
    double setpoint,
    const SteadySetpointOptions& options,
    const std::function<SteadyPlantMeasurements(double)>& evaluate_plant)
{
    if (!evaluate_plant ||
        !(std::isfinite(setpoint) && setpoint > 0.0 &&
          std::isfinite(options.minimum_speed_rpm) &&
          std::isfinite(options.maximum_speed_rpm) &&
          options.minimum_speed_rpm > 0.0 &&
          options.maximum_speed_rpm > options.minimum_speed_rpm &&
          !std::isnan(options.maximum_shaft_power_w) &&
          options.maximum_shaft_power_w > 0.0 &&
          std::isfinite(options.absolute_tolerance) &&
          options.absolute_tolerance >= 0.0 &&
          std::isfinite(options.relative_tolerance) &&
          options.relative_tolerance > 0.0 &&
          options.maximum_iterations > 0)) {
        throw std::invalid_argument(
            "Invalid steady compressor-control request");
    }

    const auto variable = [mode](const SteadyPlantMeasurements& m) {
        switch (mode) {
        case SteadySetpointMode::SuctionPressure:
            return m.suction_pressure_pa;
        case SteadySetpointMode::DischargePressure:
            return m.discharge_pressure_pa;
        case SteadySetpointMode::PipelineExportFlow:
            return m.pipeline_export_flow_kg_s;
        case SteadySetpointMode::ManualSpeed:
            break;
        }
        throw std::invalid_argument(
            "Manual speed is not a root-solved process variable");
    };

    SteadySetpointResult result;
    const bool power_cap_enabled =
        std::isfinite(options.maximum_shaft_power_w);
    auto plant_at = [&](double speed) {
        auto measurements = evaluate_plant(speed);
        ++result.plant_evaluations;
        validate_plant_measurements(measurements);
        if (power_cap_enabled &&
            !(std::isfinite(measurements.compressor_shaft_power_w) &&
              measurements.compressor_shaft_power_w >= 0.0)) {
            throw std::runtime_error(
                "Plant evaluator must report finite nonnegative shaft power "
                "when a power limit is enabled");
        }
        return measurements;
    };
    auto residual_at = [&](const SteadyPlantMeasurements& m) {
        return variable(m) - setpoint;
    };

    constexpr int power_search_iterations = 80;
    auto highest_power_feasible = [&](double feasible_speed,
                                      SteadyPlantMeasurements feasible,
                                      double infeasible_speed) {
        for (int i = 0; i < power_search_iterations; ++i) {
            const double mid = feasible_speed +
                0.5 * (infeasible_speed - feasible_speed);
            if (mid == feasible_speed || mid == infeasible_speed) break;
            auto trial = plant_at(mid);
            if (trial.compressor_shaft_power_w <=
                options.maximum_shaft_power_w) {
                feasible_speed = mid;
                feasible = std::move(trial);
            } else {
                infeasible_speed = mid;
            }
        }
        return std::pair{feasible_speed, std::move(feasible)};
    };

    double low = options.minimum_speed_rpm;
    double high = options.maximum_speed_rpm;

    if (mode == SteadySetpointMode::ManualSpeed) {
        const double requested_speed = setpoint;
        double speed = std::clamp(
            requested_speed, low, high);
        auto measurements = plant_at(speed);
        result.speed_rpm = speed;
        result.measured_value = speed;
        result.measurements = measurements;

        const bool speed_clamped = speed != requested_speed;
        if (power_cap_enabled &&
            measurements.compressor_shaft_power_w >
                options.maximum_shaft_power_w) {
            auto minimum_measurements = speed == low
                ? measurements : plant_at(low);
            if (minimum_measurements.compressor_shaft_power_w >
                options.maximum_shaft_power_w) {
                result.speed_rpm = low;
                result.measured_value = low;
                result.measurements = std::move(minimum_measurements);
                result.status = SteadySetpointStatus::PowerLimitReached;
                return result;
            }
            const auto feasible = highest_power_feasible(
                low, std::move(minimum_measurements), speed);
            result.speed_rpm = feasible.first;
            result.measured_value = feasible.first;
            result.measurements = feasible.second;
            result.status = SteadySetpointStatus::PowerLimitReached;
            return result;
        }
        result.status = speed_clamped
            ? SteadySetpointStatus::SpeedLimitReached
            : SteadySetpointStatus::Converged;
        return result;
    }

    auto low_measurements = plant_at(low);
    if (power_cap_enabled &&
        low_measurements.compressor_shaft_power_w >
            options.maximum_shaft_power_w) {
        result.status = SteadySetpointStatus::PowerLimitReached;
        result.speed_rpm = low;
        result.measured_value = variable(low_measurements);
        result.residual = result.measured_value - setpoint;
        result.measurements = std::move(low_measurements);
        return result;
    }

    auto high_measurements = plant_at(high);
    bool power_limited = false;
    if (power_cap_enabled &&
        high_measurements.compressor_shaft_power_w >
            options.maximum_shaft_power_w) {
        const auto feasible = highest_power_feasible(
            low, std::move(low_measurements), high);
        high = feasible.first;
        high_measurements = feasible.second;
        power_limited = true;
    }

    auto low_eval = std::pair{residual_at(low_measurements), low_measurements};
    auto high_eval = std::pair{residual_at(high_measurements), high_measurements};
    const double tolerance = options.absolute_tolerance +
        options.relative_tolerance * std::abs(setpoint);

    if (std::abs(low_eval.first) <= tolerance) {
        result.status = SteadySetpointStatus::Converged;
        result.speed_rpm = low;
        result.residual = low_eval.first;
        result.measured_value = low_eval.first + setpoint;
        result.measurements = low_eval.second;
        return result;
    }
    if (std::abs(high_eval.first) <= tolerance) {
        result.status = SteadySetpointStatus::Converged;
        result.speed_rpm = high;
        result.residual = high_eval.first;
        result.measured_value = high_eval.first + setpoint;
        result.measurements = high_eval.second;
        return result;
    }

    if (std::signbit(low_eval.first) == std::signbit(high_eval.first)) {
        const bool use_low =
            std::abs(low_eval.first) <= std::abs(high_eval.first);
        result.status = power_limited && !use_low
            ? SteadySetpointStatus::PowerLimitReached
            : SteadySetpointStatus::TargetNotBracketed;
        result.speed_rpm = use_low ? low : high;
        result.residual = use_low ? low_eval.first : high_eval.first;
        result.measured_value = result.residual + setpoint;
        result.measurements = use_low ? low_eval.second : high_eval.second;
        return result;
    }

    for (int iteration = 0;
         iteration < options.maximum_iterations;
         ++iteration) {
        const double mid = low + 0.5 * (high - low);
        auto measurements = plant_at(mid);
        if (power_cap_enabled &&
            measurements.compressor_shaft_power_w >
                options.maximum_shaft_power_w) {
            high = mid;
            high_eval = {residual_at(measurements), std::move(measurements)};
            power_limited = true;
            continue;
        }
        auto mid_eval = std::pair{
            residual_at(measurements), std::move(measurements)};
        result.speed_rpm = mid;
        result.residual = mid_eval.first;
        result.measured_value = mid_eval.first + setpoint;
        result.measurements = mid_eval.second;

        if (std::abs(mid_eval.first) <= tolerance) {
            result.status = SteadySetpointStatus::Converged;
            return result;
        }

        if (std::signbit(mid_eval.first) == std::signbit(low_eval.first)) {
            low = mid;
            low_eval = std::move(mid_eval);
        } else {
            high = mid;
            high_eval = std::move(mid_eval);
        }
    }

    result.status = power_limited
        ? SteadySetpointStatus::PowerLimitReached
        : SteadySetpointStatus::IterationLimit;
    return result;
}

// A3: isenthalpic throttling to suction pressure. An optional cooler model
// may be applied after this function, before the suction mixer.
inline thermomodel::State isenthalpic_recycle_return(
    const thermomodel::Model& model,
    const thermomodel::State& discharge_state,
    double suction_pressure_pa)
{
    if (!(std::isfinite(discharge_state.h) &&
          std::isfinite(discharge_state.T) &&
          std::isfinite(suction_pressure_pa) &&
          suction_pressure_pa > 0.0 &&
          suction_pressure_pa < discharge_state.p)) {
        throw std::invalid_argument(
            "Invalid discharge state or recycle suction pressure");
    }
    return invert_temperature(
        model,
        suction_pressure_pa,
        discharge_state.h,
        Property::Enthalpy,
        discharge_state.T).state;
}

inline thermomodel::State mix_suction_streams(
    const thermomodel::Model& model,
    const CompressorPipePort& fresh_inlet,
    double fresh_flow_kg_s,
    const thermomodel::State& recycle_return,
    double recycle_flow_kg_s,
    double pressure_tolerance = 1.0e-8)
{
    if (!(std::isfinite(fresh_flow_kg_s) && fresh_flow_kg_s > 0.0 &&
          std::isfinite(recycle_flow_kg_s) && recycle_flow_kg_s >= 0.0 &&
          std::isfinite(pressure_tolerance) && pressure_tolerance > 0.0)) {
        throw std::invalid_argument("Invalid suction mixing inputs");
    }
    const auto fresh = evaluate(
        model, fresh_inlet.static_temperature_k,
        fresh_inlet.static_pressure_pa);
    if (recycle_flow_kg_s == 0.0) {
        return fresh;
    }
    if (std::abs(recycle_return.p - fresh_inlet.static_pressure_pa) >
        pressure_tolerance *
            std::max(1.0, fresh_inlet.static_pressure_pa)) {
        throw std::invalid_argument(
            "Recycle return state must be at the suction mixer pressure");
    }

    const double total_flow = fresh_flow_kg_s + recycle_flow_kg_s;
    const double mixed_enthalpy =
        (fresh_flow_kg_s * fresh.h +
         recycle_flow_kg_s * recycle_return.h) / total_flow;
    return invert_temperature(
        model,
        fresh_inlet.static_pressure_pa,
        mixed_enthalpy,
        Property::Enthalpy,
        fresh.T).state;
}

struct RecycleValveDomain {
    double minimum_suction_pressure_pa = 0.0;
    double maximum_suction_pressure_pa = 0.0;
    double minimum_discharge_pressure_pa = 0.0;
    double maximum_discharge_pressure_pa = 0.0;
    double maximum_flow_kg_s = 0.0;
};

struct RecycleValveMap {
    RecycleValveDomain domain;
    std::function<double(
        double opening,
        double suction_pressure_pa,
        const thermomodel::State& discharge_state)> calibrated_flow_kg_s;

    void validate() const
    {
        const auto& d = domain;
        if (!(std::isfinite(d.minimum_suction_pressure_pa) &&
              std::isfinite(d.maximum_suction_pressure_pa) &&
              d.minimum_suction_pressure_pa > 0.0 &&
              d.maximum_suction_pressure_pa > d.minimum_suction_pressure_pa &&
              std::isfinite(d.minimum_discharge_pressure_pa) &&
              std::isfinite(d.maximum_discharge_pressure_pa) &&
              d.minimum_discharge_pressure_pa > 0.0 &&
              d.maximum_discharge_pressure_pa > d.minimum_discharge_pressure_pa &&
              std::isfinite(d.maximum_flow_kg_s) && d.maximum_flow_kg_s > 0.0)) {
            throw std::invalid_argument("Invalid recycle-valve calibration domain");
        }
    }

    double flow(
        double opening,
        double suction_pressure_pa,
        const thermomodel::State& discharge_state) const
    {
        validate();
        if (!(std::isfinite(opening) && opening >= 0.0 && opening <= 1.0 &&
              std::isfinite(suction_pressure_pa) && suction_pressure_pa > 0.0 &&
              std::isfinite(discharge_state.p) &&
              std::isfinite(discharge_state.h) &&
              std::isfinite(discharge_state.T))) {
            throw std::invalid_argument("Invalid recycle-valve query");
        }
        if (opening == 0.0) return 0.0;
        if (!calibrated_flow_kg_s ||
            suction_pressure_pa < domain.minimum_suction_pressure_pa ||
            suction_pressure_pa > domain.maximum_suction_pressure_pa ||
            discharge_state.p < domain.minimum_discharge_pressure_pa ||
            discharge_state.p > domain.maximum_discharge_pressure_pa ||
            !(discharge_state.p > suction_pressure_pa)) {
            throw std::out_of_range(
                "Recycle-valve query is outside its calibrated domain");
        }
        const double value = calibrated_flow_kg_s(
            opening, suction_pressure_pa, discharge_state);
        if (!(std::isfinite(value) && value >= 0.0 &&
              value <= domain.maximum_flow_kg_s)) {
            throw std::runtime_error(
                "Recycle valve map returned an invalid flow");
        }
        return value;
    }
};

struct RecycleSteadyOptions {
    int maximum_iterations = 60;
    double relaxation = 0.5;
    double absolute_flow_tolerance_kg_s = 1.0e-7;
    double relative_flow_tolerance = 1.0e-6;
    double relative_return_enthalpy_tolerance = 1.0e-7;
};

struct RecycleSteadyResult {
    bool converged = false;
    std::string message;
    int iterations = 0;
    double pipeline_export_flow_kg_s = 0.0;
    double compressor_throughflow_kg_s = 0.0;
    double recycle_flow_kg_s = 0.0;
    double recycle_mass_residual_kg_s = 0.0;
    double recycle_cooler_heat_removed_w = 0.0;
    double suction_mixer_enthalpy_residual_w = 0.0;
    CompressorPipePort mixed_suction;
    CompressorComponentResult compressor;
};

template<class CompressorEvaluator, class RecycleReturnModel>
RecycleSteadyResult solve_steady_recycle(
    const thermomodel::Model& model,
    const CompressorPipePort& fresh_inlet,
    double pipeline_export_flow_kg_s,
    double valve_opening,
    const RecycleValveMap& valve_map,
    CompressorEvaluator&& evaluate_compressor,
    RecycleReturnModel&& return_model,
    const RecycleSteadyOptions& options = {})
{
    if (!(std::isfinite(pipeline_export_flow_kg_s) &&
          pipeline_export_flow_kg_s > 0.0 &&
          std::isfinite(valve_opening) &&
          valve_opening >= 0.0 && valve_opening <= 1.0 &&
          options.maximum_iterations > 0 &&
          std::isfinite(options.relaxation) &&
          options.relaxation > 0.0 && options.relaxation <= 1.0 &&
          std::isfinite(options.absolute_flow_tolerance_kg_s) &&
          options.absolute_flow_tolerance_kg_s > 0.0 &&
          std::isfinite(options.relative_flow_tolerance) &&
          options.relative_flow_tolerance > 0.0 &&
          std::isfinite(options.relative_return_enthalpy_tolerance) &&
          options.relative_return_enthalpy_tolerance > 0.0)) {
        throw std::invalid_argument("Invalid steady recycle request");
    }

    RecycleSteadyResult result;
    result.pipeline_export_flow_kg_s = pipeline_export_flow_kg_s;
    double recycle_flow = 0.0;
    auto recycle_state = evaluate(
        model, fresh_inlet.static_temperature_k,
        fresh_inlet.static_pressure_pa);

    for (int iteration = 1;
         iteration <= options.maximum_iterations;
         ++iteration) {
        const auto mixed = mix_suction_streams(
            model, fresh_inlet, pipeline_export_flow_kg_s,
            recycle_state, recycle_flow);
        const double throughflow = pipeline_export_flow_kg_s + recycle_flow;
        auto compressor = evaluate_compressor(mixed, throughflow);
        result.iterations = iteration;
        result.compressor = std::move(compressor);
        result.compressor_throughflow_kg_s = throughflow;
        result.recycle_flow_kg_s = recycle_flow;
        result.mixed_suction = {
            mixed.T, fresh_inlet.static_pressure_pa
        };

        if (!result.compressor.valid ||
            !result.compressor.outlet_static_state) {
            result.message = "Compressor evaluation failed in recycle loop";
            return result;
        }

        if (valve_opening == 0.0) {
            result.converged = true;
            result.message = "Recycle valve closed";
            return result;
        }

        const auto& discharge = *result.compressor.outlet_static_state;
        const double requested_flow = valve_map.flow(
            valve_opening, fresh_inlet.static_pressure_pa, discharge);
        auto requested_return = return_model(
            discharge, fresh_inlet.static_pressure_pa);
        if (!(std::isfinite(requested_return.T) && requested_return.T > 0.0 &&
              std::isfinite(requested_return.p) && requested_return.p > 0.0 &&
              std::isfinite(requested_return.h) &&
              std::abs(requested_return.p - fresh_inlet.static_pressure_pa) <=
                  1.0e-8 * std::max(1.0, fresh_inlet.static_pressure_pa))) {
            throw std::runtime_error(
                "Recycle return model must produce a finite suction-pressure state");
        }
        const auto checked_return = turbocompgpu::evaluate(
            model, requested_return.T, requested_return.p);
        if (std::abs(checked_return.h - requested_return.h) >
            1.0e-8 * std::max({1.0, std::abs(checked_return.h),
                               std::abs(requested_return.h)})) {
            throw std::runtime_error(
                "Recycle return state is inconsistent with the selected EOS");
        }
        requested_return = checked_return;
        const double flow_tolerance =
            options.absolute_flow_tolerance_kg_s +
            options.relative_flow_tolerance *
                std::max(requested_flow, recycle_flow);
        const double h_scale =
            std::max({1.0, std::abs(requested_return.h),
                      std::abs(recycle_state.h)});
        const double h_error =
            std::abs(requested_return.h - recycle_state.h);

        result.recycle_mass_residual_kg_s =
            requested_flow - recycle_flow;
        result.recycle_cooler_heat_removed_w =
            requested_flow * (discharge.h - requested_return.h);

        const double mixed_flow = pipeline_export_flow_kg_s + recycle_flow;
        const auto checked_mix = mix_suction_streams(
            model, fresh_inlet, pipeline_export_flow_kg_s,
            recycle_state, recycle_flow);
        result.suction_mixer_enthalpy_residual_w =
            mixed_flow * checked_mix.h -
            pipeline_export_flow_kg_s *
                evaluate(model, fresh_inlet.static_temperature_k,
                         fresh_inlet.static_pressure_pa).h -
            recycle_flow * recycle_state.h;

        if (std::abs(result.recycle_mass_residual_kg_s) <= flow_tolerance &&
            h_error <= options.relative_return_enthalpy_tolerance * h_scale) {
            result.converged = true;
            result.message = "Recycle mass/return-state iteration converged";
            return result;
        }

        recycle_flow += options.relaxation *
            result.recycle_mass_residual_kg_s;
        recycle_flow = std::clamp(
            recycle_flow, 0.0, valve_map.domain.maximum_flow_kg_s);
        recycle_state = std::move(requested_return);
    }

    result.message = "Recycle fixed-point iteration did not converge";
    return result;
}

} // namespace turbocompgpu
