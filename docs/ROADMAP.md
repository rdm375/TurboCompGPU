# TurboCompGPU Development Roadmap

Reviewed: 2026-11-11
Code status reference: branch main, base commit bcb7b89 (2026-11-06)

## Project goals

The project has two intended outputs:

1. **Pipeline component:** a compressor model that can be connected to a transient pipe solver and participate in pipeline simulation.
2. **Semi-realistic surrogate candidate:** a fast ReLU interpolant of sampled compressor operation, continuing the approach explored in the original NN Interpolation hackathon.

These goals share data and interfaces, but they are not the same deliverable. The compressor physics solver is the reference model and possible sample generator. The surrogate is a fitted approximation with its own accuracy, domain, and validation requirements.

## Current status

| Area | Status | Evidence and limitation |
|---|---|---|
| Thermodynamics and nonlinear compressor solver | Implemented and numerically tested | IdealGas and PengRobinson paths, conservation checks, reference-invariance and solver validation |
| Operating points and maps | Implemented for tested domains | Deterministic sweeps and continuation; tested convergence is not proof of physical admissibility |
| Impeller/diffuser complete-stage model | Implemented with numerical verification | Frictional and variable-width diffuser work is documented; empirical calibration is still missing |
| Incidence closure | Experimental model extension | 105 paired points across three thermodynamic configurations; coefficient not experimentally calibrated |
| Diffuser profile geometry | Infrastructure present | PCHIP hub/shroud profiles and tests; actual NASA HECC hot-wall geometry is not reconstructed |
| Pipeline coupling | **Not implemented** | No end-to-end compressor boundary connected to a transient pipe model |
| ReLU surrogate | **Preliminary baseline added; not yet qualified** | A trainer fits current C++-generated samples; its structured holdout results have not yet been recorded, and it is not connected to the pipeline |
| GPU execution | **Not implemented** | Current code and reported timings are CPU based |
| Experimental validation | **Air screening only; natural-gas validation not complete** | Eckardt A/B and Choi–Kang provide limited air-compressor evidence. The current Eckardt runner supports ideal gas and PR air, not a specified natural-gas mixture. Exit criteria and the accelerator reference gate are in [`PHYSICS_VALIDATION_EXIT_GATE.md`](PHYSICS_VALIDATION_EXIT_GATE.md). |

For supporting details, see [Project Status](PROJECT_STATUS.md) and the linked milestone reports.

## Workstream A — Connect a compressor to a pipeline

### A0 — Define the component contract

Use the existing flow-station and operating-point concepts to define a small, explicit component interface.

The contract should identify:

- Inlet state and flow variables supplied by the pipe model.
- Compressor controls and configuration, including rotational speed.
- Outlet state and shaft-work or power outputs returned to the pipe/station model.
- Solver status, diagnostics, and out-of-domain or invalid-state behavior.
- The thermodynamic property model and static/total-state conventions at both boundaries.

**Acceptance:** one documented API or adapter with unit tests for state conversion, mass-flow sign/convention, and energy consistency. Keep the interface independent of a GPU backend and avoid requiring a full station-control system in this first step.

### A1 — Steady pipe–compressor–pipe case

Place one compressor between upstream and downstream pipe sections. First establish a reproducible steady operating point and verify that the coupled pressure and flow agree with the component’s standalone result.

**Acceptance:** a small deterministic example with convergence diagnostics and mass/energy balance checks at the component boundary.

### A2 — Transient single-compressor case

Connect the same component to the transient pipe solver. Exercise a controlled change in demand or compressor speed and verify time-step behavior, coupling convergence, and conservation.

**Acceptance:** a documented transient case with repeatable results, explicit failure handling, and comparison against a tighter-coupling or refined-step reference.

### A3 — Station behavior only when needed

Add parallel units, recycle/bypass, controls, and station-level operating logic only when the pipeline use case requires them. Keep those features distinct from the first single-compressor coupling.

## Workstream B — Reproduce and qualify the ReLU interpolation idea

The original NN Interpolation hackathon used ReLU networks as interpolation functions over sampled compressor-operation input/output data. It treated the compressor simulator as a source of samples rather than requiring access to the simulator internals. Preserve that as the baseline method and record the provenance and units of every sample set.

### B0 — Establish a model-generated ReLU baseline

The original hackathon artifacts are not currently available. Use the existing compressor-map samples to recreate the core ReLU interpolation method, while labeling the data as outputs of the current C++ model rather than legacy-simulator or experimental truth. The initial script is documented in RELU_SURROGATE_BASELINE.md.

**Acceptance:** run the training script, save its structured speed-line holdout metrics, and confirm that the model checkpoints can reproduce predictions from the report.

### B1 — Define a semi-realistic training domain

Use the physics-based compressor model to generate additional samples only within the model’s documented validity domain. Where legacy simulator samples exist, compare the two sources at matched conditions instead of silently treating either as experimental truth.

Vary relevant operating inputs (for example inlet state, mass flow, and speed) while preserving geometry and gas composition metadata.

**Acceptance:** a dataset manifest, domain coverage report, solver-failure counts, and a clear label for each source: legacy simulator, current physics model, or experiment.

### B2 — Validate the interpolant

Compare held-out predictions against their sample source and against independent physics-solver evaluations. Report maximum, percentile, and RMS errors per output, boundary behavior, and out-of-domain behavior. Check physical consistency, including positive absolute pressure and nonnegative power where applicable.

**Acceptance:** a qualified domain and explicit error gate chosen for the pipeline use case. Do not use random train/test splits alone when nearby operating points can leak across the split; hold out structured regions or complete speed lines as well.

### B3 — Use the surrogate in a pipeline case

Expose the surrogate through the same component contract as the reference solver. Compare coupled pipeline outputs against the reference-model coupling over the qualified domain.

**Acceptance:** repeatable pipeline-level comparisons and a defined response when queries leave the surrogate’s qualified domain.

## Model credibility work

These items improve the reference model but do not block the first software integration prototype:

- Download, inventory, and process the NASA HECC public archive; record its checksum and use its independent multi-speed data for the campaign in [`VALIDATION_CAMPAIGN.md`](VALIDATION_CAMPAIGN.md). Earlier 2026-10-13 access attempts returned HTTP 403; the current official archive URL and guide have since been located, but local retrieval and extraction still need confirmation.
- Reconstruct actual hub/shroud geometry and compare diffuser predictions with experimental data.
- Replace exploratory or prescribed loss closures with calibrated or independently supported models.
- Establish operating validity limits. The current work does not provide a predictive surge boundary; distinguish solver convergence from physical admissibility.
- Resolve the intentionally disabled test and record a fresh full regression result for the latest source tree.

## Acceleration

Keep a reproducible scalar CPU build as the oracle for backend parity. Accelerator parity and throughput experiments may proceed against that frozen oracle, but do not claim natural-gas predictive accuracy until the natural-gas validation gate in [`PHYSICS_VALIDATION_EXIT_GATE.md`](PHYSICS_VALIDATION_EXIT_GATE.md) passes. No GPU solver is currently implemented, and past CPU timings must not be presented as GPU performance.

## Working sequence

1. A0: define and test the component contract.
2. B0: run and review the C++-sample ReLU baseline; recover legacy hackathon artifacts later if they become available.
3. A1: run one steady coupled pipeline case.
4. B1–B2: build the semi-realistic data set and qualify the ReLU interpolant.
5. A2 and B3: run transient comparisons with both reference model and surrogate.
6. Pursue experimental validation and acceleration against identified use cases.

### Experimental validation checkpoint — 2026-11-08

The radial-discharge Eckardt Impeller O traverse has been identified as
cross-family data and is not used to qualify the backswept Impeller A model.
The new same-family screen digitizes three 14,000 rpm points from Choi and
Kang (1999), Figure 2, and compares predicted total pressure at `r/r2=1.075`
and total temperature at `r/r2=1.687` using the source's mixed-plane
efficiency definition. It gives 2.43% pressure-ratio MAPE, mean absolute
temperature-rise-coefficient error 0.0531, and 7.65 percentage-points
efficiency MAE. Treat these as a narrow screening result: the three source
symbols were digitized from a plot, and they cover only one speed line. See
[`ECKARDT_A_CHOI_KANG_1999_QUALIFICATION.md`](ECKARDT_A_CHOI_KANG_1999_QUALIFICATION.md)
and `benchmarks/results-eckardt-impeller-a/choi-kang-1999-fig2/`.
