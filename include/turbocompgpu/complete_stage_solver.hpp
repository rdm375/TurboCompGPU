#pragma once

#include <turbocompgpu/complete_stage_branches.hpp>

#include <optional>
#include <string>
#include <exception>

namespace turbocompgpu {

// Stage 8B.5a: one-call complete-stage orchestration.
//
// This layer introduces no new physical equations.
// It composes the independently validated impeller solver
// and the branch-aware lossless vaneless diffuser solver.
//
// Complete means numerical/physical consistency with the
// implemented model. It does not establish surge margin,
// compressor-system choking, or experimental accuracy.

enum class CompleteStageStatus {
    ImpellerFailed,
    DiffuserFailed,
    Complete
};

struct CompleteStageResult {
    CompleteStageStatus status =
        CompleteStageStatus::ImpellerFailed;

    DiffuserBranch branch =
        DiffuserBranch::Subsonic;

    // Always retained, including unsuccessful impeller solves.
    OperatingPointResult impeller;

    // Present once a validated impeller solution has been
    // passed to the complete-stage integration layer.
    std::optional<CompleteStageBranchResult> stage;

    std::string message;

    bool complete() const noexcept
    {
        return status == CompleteStageStatus::Complete;
    }

    bool impeller_valid() const noexcept
    {
        return impeller.converged &&
               impeller.independently_validated &&
               impeller.basic_physics_valid;
    }

    bool station3_valid() const noexcept
    {
        return stage.has_value() &&
               stage->station3.has_value() &&
               (stage->status == DiffuserBranchStatus::Converged ||
                stage->status == DiffuserBranchStatus::Critical);
    }
};

inline CompleteStageResult evaluate_complete_stage(
    const thermomodel::Model& model,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& impeller_geometry,
    const VanelessDiffuserGeometry& diffuser_geometry,
    const OperatingPointRequest& request,
    const NonlinearResidualScales& scales,
    const NonlinearSolverOptions& options = {},
    DiffuserBranch branch = DiffuserBranch::Subsonic,
    double validation_tolerance = 1.0e-8)
{
    CompleteStageResult result;
    result.branch = branch;

    // First solve and independently validate the impeller.
    result.impeller = evaluate_operating_point(
        model,
        impeller_geometry,
        request,
        scales,
        options,
        validation_tolerance);

    if (!result.impeller_valid()) {
        result.status = CompleteStageStatus::ImpellerFailed;
        result.message =
            "Impeller solve or independent validation failed: " +
            result.impeller.message;
        return result;
    }

    // Only validated impeller states enter the diffuser.
    //
    // Keep the impeller result even if diffuser evaluation
    // fails or throws because of invalid geometry or a
    // thermodynamic-domain problem.
    try {
        result.stage = evaluate_complete_stage_branch(
            model,
            derivatives,
            impeller_geometry,
            request,
            result.impeller,
            diffuser_geometry,
            branch);
    }
    catch (const std::exception& error) {
        result.status = CompleteStageStatus::DiffuserFailed;
        result.message =
            std::string("Diffuser evaluation failed: ") +
            error.what();
        return result;
    }

    if (!result.station3_valid()) {
        result.status = CompleteStageStatus::DiffuserFailed;
        result.message =
            "Diffuser did not produce a valid station 3; "
            "inspect stage.diffuser.status";
        return result;
    }

    result.status = CompleteStageStatus::Complete;
    result.message =
        "Complete-stage solution converged and validated";

    return result;
}

} // namespace turbocompgpu
