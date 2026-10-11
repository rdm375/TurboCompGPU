#pragma once

#include <turbocompgpu/thermodynamic_solver.hpp>

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

struct OperatingPoint {
    double inlet_T0;
    double inlet_p0;
    double outlet_p0;
    double mass_flow;
    double isentropic_efficiency;
};

struct CompressorResult {
    thermomodel::State inlet;
    thermomodel::State outlet_isentropic;
    thermomodel::State outlet;

    double specific_work;
    double power;

    InversionResult entropy_inversion;
    InversionResult enthalpy_inversion;
};

inline CompressorResult solve_compressor(
    const thermomodel::Model& model,
    const OperatingPoint& input)
{
    if (!std::isfinite(input.inlet_T0) ||
        !std::isfinite(input.inlet_p0) ||
        !std::isfinite(input.outlet_p0) ||
        !std::isfinite(input.mass_flow) ||
        !std::isfinite(input.isentropic_efficiency) ||
        input.inlet_T0 <= 0.0 ||
        input.inlet_p0 <= 0.0 ||
        input.outlet_p0 <= input.inlet_p0 ||
        input.mass_flow <= 0.0 ||
        input.isentropic_efficiency <= 0.0 ||
        input.isentropic_efficiency > 1.0) {
        throw std::invalid_argument(
            "Invalid compressor operating point");
    }

    CompressorResult result{};

    result.inlet = evaluate(
        model, input.inlet_T0, input.inlet_p0);

    result.entropy_inversion = invert_temperature(
        model,
        input.outlet_p0,
        result.inlet.s,
        Property::Entropy,
        input.inlet_T0);

    result.outlet_isentropic =
        result.entropy_inversion.state;

    const double delta_h_isentropic =
        result.outlet_isentropic.h - result.inlet.h;

    if (!std::isfinite(delta_h_isentropic) ||
        delta_h_isentropic <= 0.0) {
        throw std::runtime_error(
            "Non-positive isentropic enthalpy rise");
    }

    result.specific_work =
        delta_h_isentropic /
        input.isentropic_efficiency;

    const double target_h =
        result.inlet.h + result.specific_work;

    result.enthalpy_inversion = invert_temperature(
        model,
        input.outlet_p0,
        target_h,
        Property::Enthalpy,
        result.outlet_isentropic.T);

    result.outlet =
        result.enthalpy_inversion.state;

    result.power =
        input.mass_flow * result.specific_work;

    if (!std::isfinite(result.specific_work) ||
        !std::isfinite(result.power)) {
        throw std::runtime_error(
            "Non-finite compressor work or power");
    }

    return result;
}

} // namespace turbocompgpu
