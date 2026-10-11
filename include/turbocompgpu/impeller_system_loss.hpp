#pragma once

#include <turbocompgpu/impeller_system.hpp>
#include <turbocompgpu/loss.hpp>

namespace turbocompgpu {

struct ImpellerSystemWithLossResult {
    ImpellerSystemResult system;
    EntropyLossResult loss;
};

inline ImpellerSystemWithLossResult
evaluate_impeller_system_with_loss(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& geometry,
    const ImpellerTrialState& trial,
    const EntropyLossModel& loss_model)
{
    // Evaluate the complete impeller trial state.
    //
    // This calculation remains independent of the
    // prescribed entropy-generation closure.

    const auto system =
        evaluate_impeller_system(
            model,
            geometry,
            trial);

    // Evaluate entropy generation using the actual
    // inlet and outlet stagnation states.

    const auto loss =
        evaluate_entropy_loss(
            system.inlet_stagnation.stagnation_state,
            system.outlet_stagnation.stagnation_state,
            loss_model);

    return {
        system,
        loss
    };
}

} // namespace turbocompgpu
