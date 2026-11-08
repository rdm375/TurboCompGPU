#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

namespace turbocompgpu {

enum class CompressorControlMode {
    DischargePressure,
    SuctionPressure,
    PipelineMassFlow,
    ManualSpeed
};

enum class CompressorProcessVariable : std::size_t {
    DischargePressure = 0,
    SuctionPressure = 1,
    PipelineMassFlow = 2
};

// Gains are in rpm per process-variable unit and rpm per
// process-variable-unit second. direction is +1 or -1 and is
// configured for the measured plant response.
struct CompressorPiGains {
    double proportional = 0.0;
    double integral = 0.0;
    double direction = 1.0;
};

struct AntiSurgeControlConfig {
    bool enabled = false;
    double target_margin_fraction = 0.10;
    double proportional = 0.0;
    double integral = 0.0;
    double valve_rate_limit_per_second = 1.0;
};

struct CompressorControlConfig {
    double minimum_speed_rpm = 0.0;
    double maximum_speed_rpm = 0.0;
    double speed_rate_limit_rpm_per_second = 0.0;

    // Pressure units: Pa. Flow units: kg/s.
    // Tune each loop against its compressor + pipeline response.
    std::array<CompressorPiGains, 3> process_gains{};

    AntiSurgeControlConfig anti_surge{};
};

struct CompressorControlMeasurement {
    double suction_pressure_pa = 0.0;
    double discharge_pressure_pa = 0.0;

    // Net flow exported to the pipeline. With recycle active this
    // differs from the compressor's internal through-flow.
    double pipeline_mass_flow_kg_s = 0.0;

    double speed_rpm = 0.0;
    double recycle_valve_opening = 0.0;

    // Must be supplied by a separately qualified compressor surge map.
    // Positive margin is distance from the surge line, normalized by
    // surge-line flow; negative values indicate crossing that line.
    std::optional<double> surge_margin_fraction;
};

struct CompressorControlSetpoint {
    CompressorControlMode mode = CompressorControlMode::DischargePressure;
    double value = 0.0; // Pa, kg/s, or rpm according to mode.
};

struct CompressorControlCommand {
    double speed_command_rpm = 0.0;
    double recycle_valve_opening = 0.0;
    CompressorControlMode mode = CompressorControlMode::DischargePressure;
    bool speed_rate_limited = false;
    bool speed_at_minimum = false;
    bool speed_at_maximum = false;
    bool anti_surge_override_active = false;
};

class CompressorController {
public:
    explicit CompressorController(CompressorControlConfig config)
        : config_(std::move(config))
    {
        validate_config();
    }

    CompressorControlCommand step(
        const CompressorControlSetpoint& setpoint,
        const CompressorControlMeasurement& measurement,
        double time_step_seconds)
    {
        validate_measurement(measurement);
        if (!(std::isfinite(time_step_seconds) &&
              time_step_seconds > 0.0)) {
            throw std::invalid_argument(
                "Controller time step must be positive and finite");
        }
        if (!(std::isfinite(setpoint.value) && setpoint.value > 0.0)) {
            throw std::invalid_argument(
                "Compressor setpoint must be positive and finite");
        }
        if (config_.anti_surge.enabled &&
            (!measurement.surge_margin_fraction.has_value() ||
             !std::isfinite(*measurement.surge_margin_fraction))) {
            throw std::invalid_argument(
                "Anti-surge control requires a qualified surge-margin input");
        }

        const double prior_speed = initialized_
            ? last_speed_command_rpm_
            : clamp(measurement.speed_rpm,
                    config_.minimum_speed_rpm,
                    config_.maximum_speed_rpm);

        double speed_command = prior_speed;
        bool speed_rate_limited = false;

        if (setpoint.mode == CompressorControlMode::ManualSpeed) {
            speed_command = limit_speed(
                setpoint.value, prior_speed, time_step_seconds,
                speed_rate_limited);
            last_speed_error_ = 0.0;
            last_mode_ = setpoint.mode;
        } else {
            const auto variable = process_variable(setpoint.mode);
            const auto index = static_cast<std::size_t>(variable);
            const auto& gains = config_.process_gains[index];

            if (gains.proportional == 0.0 && gains.integral == 0.0) {
                throw std::invalid_argument(
                    "Selected compressor control loop has zero gains");
            }

            const double process_value =
                measured_value(variable, measurement);
            const double error = setpoint.value - process_value;
            const double previous_error =
                (!initialized_ || last_mode_ != setpoint.mode)
                    ? 0.0
                    : last_speed_error_;

            const double proportional_increment =
                gains.proportional * (error - previous_error);
            const double integral_increment =
                gains.integral * error * time_step_seconds;

            double raw = prior_speed + gains.direction *
                (proportional_increment + integral_increment);

            speed_command = limit_speed(
                raw, prior_speed, time_step_seconds,
                speed_rate_limited);

            // Conditional integration prevents the integral term from
            // winding farther into either a speed or slew-rate limit.
            const bool integrates_into_limit =
                (raw > speed_command && gains.direction *
                    integral_increment > 0.0) ||
                (raw < speed_command && gains.direction *
                    integral_increment < 0.0);
            if (integrates_into_limit) {
                raw = prior_speed +
                    gains.direction * proportional_increment;
                speed_command = limit_speed(
                    raw, prior_speed, time_step_seconds,
                    speed_rate_limited);
            }

            last_speed_error_ = error;
            last_mode_ = setpoint.mode;
        }

        const double recycle_command = update_recycle(
            measurement, time_step_seconds);

        initialized_ = true;
        last_speed_command_rpm_ = speed_command;
        last_recycle_command_ = recycle_command;

        return {
            speed_command,
            recycle_command,
            setpoint.mode,
            speed_rate_limited,
            speed_command <= config_.minimum_speed_rpm,
            speed_command >= config_.maximum_speed_rpm,
            config_.anti_surge.enabled &&
                measurement.surge_margin_fraction.has_value() &&
                *measurement.surge_margin_fraction <=
                    config_.anti_surge.target_margin_fraction
        };
    }

    void reset() noexcept
    {
        initialized_ = false;
        last_speed_command_rpm_ = 0.0;
        last_recycle_command_ = 0.0;
        last_speed_error_ = 0.0;
        last_surge_error_ = 0.0;
        last_mode_ = CompressorControlMode::DischargePressure;
    }

private:
    static double clamp(double value, double low, double high)
    {
        return std::max(low, std::min(value, high));
    }

    static CompressorProcessVariable process_variable(
        CompressorControlMode mode)
    {
        switch (mode) {
        case CompressorControlMode::DischargePressure:
            return CompressorProcessVariable::DischargePressure;
        case CompressorControlMode::SuctionPressure:
            return CompressorProcessVariable::SuctionPressure;
        case CompressorControlMode::PipelineMassFlow:
            return CompressorProcessVariable::PipelineMassFlow;
        case CompressorControlMode::ManualSpeed:
            break;
        }
        throw std::invalid_argument(
            "Manual speed has no automatic process variable");
    }

    static double measured_value(
        CompressorProcessVariable variable,
        const CompressorControlMeasurement& measurement)
    {
        switch (variable) {
        case CompressorProcessVariable::DischargePressure:
            return measurement.discharge_pressure_pa;
        case CompressorProcessVariable::SuctionPressure:
            return measurement.suction_pressure_pa;
        case CompressorProcessVariable::PipelineMassFlow:
            return measurement.pipeline_mass_flow_kg_s;
        }
        throw std::invalid_argument("Unknown compressor process variable");
    }

    double limit_speed(
        double requested,
        double previous,
        double time_step_seconds,
        bool& rate_limited) const
    {
        const double step = config_.speed_rate_limit_rpm_per_second *
            time_step_seconds;
        const double rate_low = std::max(
            config_.minimum_speed_rpm, previous - step);
        const double rate_high = std::min(
            config_.maximum_speed_rpm, previous + step);
        const double result = clamp(requested, rate_low, rate_high);
        rate_limited = result != clamp(
            requested,
            config_.minimum_speed_rpm,
            config_.maximum_speed_rpm);
        return result;
    }

    double update_recycle(
        const CompressorControlMeasurement& measurement,
        double time_step_seconds)
    {
        const double previous = initialized_
            ? last_recycle_command_
            : measurement.recycle_valve_opening;

        if (!config_.anti_surge.enabled) {
            const double step =
                config_.anti_surge.valve_rate_limit_per_second *
                time_step_seconds;
            return clamp(0.0, std::max(0.0, previous - step),
                         std::min(1.0, previous + step));
        }

        if (!measurement.surge_margin_fraction.has_value() ||
            !std::isfinite(*measurement.surge_margin_fraction)) {
            throw std::invalid_argument(
                "Anti-surge control requires a qualified surge-margin input");
        }

        const double error = config_.anti_surge.target_margin_fraction -
            *measurement.surge_margin_fraction;
        const double previous_error =
            initialized_ ? last_surge_error_ : 0.0;
        const double raw = previous +
            config_.anti_surge.proportional *
                (error - previous_error) +
            config_.anti_surge.integral * error * time_step_seconds;

        const double step =
            config_.anti_surge.valve_rate_limit_per_second *
            time_step_seconds;
        const double low = std::max(0.0, previous - step);
        const double high = std::min(1.0, previous + step);
        double command = clamp(raw, low, high);

        const bool integral_pushes_further =
            (raw > command &&
             config_.anti_surge.integral * error > 0.0) ||
            (raw < command &&
             config_.anti_surge.integral * error < 0.0);
        if (integral_pushes_further) {
            const double without_integral = previous +
                config_.anti_surge.proportional *
                    (error - previous_error);
            command = clamp(without_integral, low, high);
        }

        last_surge_error_ = error;
        return command;
    }

    void validate_config() const
    {
        if (!(std::isfinite(config_.minimum_speed_rpm) &&
              std::isfinite(config_.maximum_speed_rpm) &&
              config_.minimum_speed_rpm > 0.0 &&
              config_.maximum_speed_rpm >
                  config_.minimum_speed_rpm &&
              std::isfinite(
                  config_.speed_rate_limit_rpm_per_second) &&
              config_.speed_rate_limit_rpm_per_second > 0.0)) {
            throw std::invalid_argument(
                "Invalid compressor speed or ramp limits");
        }

        for (const auto& gains : config_.process_gains) {
            if (!(std::isfinite(gains.proportional) &&
                  gains.proportional >= 0.0 &&
                  std::isfinite(gains.integral) &&
                  gains.integral >= 0.0 &&
                  (gains.direction == 1.0 ||
                   gains.direction == -1.0))) {
                throw std::invalid_argument(
                    "Invalid compressor PI gains or direction");
            }
        }

        const auto& surge = config_.anti_surge;
        if (!(std::isfinite(surge.target_margin_fraction) &&
              surge.target_margin_fraction >= 0.0 &&
              std::isfinite(surge.proportional) &&
              surge.proportional >= 0.0 &&
              std::isfinite(surge.integral) &&
              surge.integral >= 0.0 &&
              std::isfinite(surge.valve_rate_limit_per_second) &&
              surge.valve_rate_limit_per_second > 0.0 &&
              (!surge.enabled ||
               surge.proportional > 0.0 ||
               surge.integral > 0.0))) {
            throw std::invalid_argument(
                "Invalid anti-surge controller configuration");
        }
    }

    static void validate_measurement(
        const CompressorControlMeasurement& measurement)
    {
        if (!(std::isfinite(measurement.suction_pressure_pa) &&
              measurement.suction_pressure_pa > 0.0 &&
              std::isfinite(measurement.discharge_pressure_pa) &&
              measurement.discharge_pressure_pa > 0.0 &&
              std::isfinite(measurement.pipeline_mass_flow_kg_s) &&
              measurement.pipeline_mass_flow_kg_s >= 0.0 &&
              std::isfinite(measurement.speed_rpm) &&
              measurement.speed_rpm > 0.0 &&
              std::isfinite(measurement.recycle_valve_opening) &&
              measurement.recycle_valve_opening >= 0.0 &&
              measurement.recycle_valve_opening <= 1.0)) {
            throw std::invalid_argument(
                "Invalid compressor control measurement");
        }
    }

    CompressorControlConfig config_;
    bool initialized_ = false;
    double last_speed_command_rpm_ = 0.0;
    double last_recycle_command_ = 0.0;
    double last_speed_error_ = 0.0;
    double last_surge_error_ = 0.0;
    CompressorControlMode last_mode_ =
        CompressorControlMode::DischargePressure;
};

} // namespace turbocompgpu
