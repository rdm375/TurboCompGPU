#pragma once

#include <thermomodel/model.hpp>

#include <cmath>
#include <stdexcept>

namespace turbocompgpu {

// Thermodynamic reference-convention adapter.
//
// Adds constants to specific enthalpy, internal energy,
// and entropy without changing physical state properties.
//
// The underlying model must outlive this adapter.

class ReferenceShiftModel final : public thermomodel::Model {
public:
    ReferenceShiftModel(
        const thermomodel::Model& base,
        double enthalpy_shift,
        double entropy_shift)
        : base_(base),
          enthalpy_shift_(enthalpy_shift),
          entropy_shift_(entropy_shift)
    {
        if (!std::isfinite(enthalpy_shift_) ||
            !std::isfinite(entropy_shift_)) {
            throw std::invalid_argument(
                "Reference shifts must be finite");
        }
    }

    thermomodel::Result<thermomodel::State>
    evaluate_Tp(
        double T,
        double pressure,
        const thermomodel::Composition& composition = {}
    ) const override
    {
        return apply_reference_shift(
            base_.evaluate_Tp(T, pressure, composition));
    }

    thermomodel::Result<thermomodel::State>
    evaluate_Trho(
        double T,
        double rho,
        const thermomodel::Composition& composition = {}
    ) const override
    {
        return apply_reference_shift(
            base_.evaluate_Trho(T, rho, composition));
    }

    thermomodel::Capabilities
    capabilities() const noexcept override
    {
        return base_.capabilities();
    }

    thermomodel::Metadata metadata() const override
    {
        auto result = base_.metadata();

        result.reference_convention +=
            " + constant enthalpy/entropy reference shifts";

        return result;
    }

private:
    thermomodel::Result<thermomodel::State>
    apply_reference_shift(
        thermomodel::Result<thermomodel::State> result
    ) const
    {
        if (!result) {
            return result;
        }

        result.value.h += enthalpy_shift_;
        result.value.u += enthalpy_shift_;
        result.value.s += entropy_shift_;

        return result;
    }

    const thermomodel::Model& base_;

    double enthalpy_shift_;
    double entropy_shift_;
};

} // namespace turbocompgpu
