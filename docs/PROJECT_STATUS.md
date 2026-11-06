# TurboCompGPU Project Status

Reviewed: 2026-11-06
Repository: [rdm375/TurboCompGPU](https://github.com/rdm375/TurboCompGPU), branch main
Latest commit observed: acaa714, “Add verified PCHIP diffuser wall geometry and HECC handoff” (2026-10-13)

## Intended outcome

TurboCompGPU is meant to support both:

1. A compressor component connected to a transient pipe model for pipeline simulation.
2. A semi-realistic surrogate candidate that continues the original NN Interpolation hackathon’s ReLU interpolation of sampled compressor operation.

The hackathon’s key idea was to learn a function from simulator input/output samples, without needing to reproduce the simulator internals. The current C++ project is a physics-based compressor model, not a completed surrogate. It can become a source of structured samples and an independent comparison model, but its own empirical validity remains limited.

## Where the code is

### Implemented and numerically exercised

- Ideal-gas and Peng–Robinson thermodynamic paths using ThermoGPU.
- Meanline impeller aerodynamics and nonlinear operating-point solves.
- Flow-station state and velocity representations, including stagnation-state evaluation.
- Deterministic compressor operating sweeps, map quantities, and continuation.
- Impeller plus vaneless-diffuser complete-stage calculations, including frictional and variable-width formulations.
- PCHIP hub/shroud profile infrastructure for variable diffuser geometry.
- Incidence-loss experiments with independent algebraic and zero-incidence checks.

The M1.5 report documents 140 converged and independently validated solves over 70 tested conditions. The M1.6b.2 report documents 105 paired operating points across three thermodynamic configurations and 315 successful nonlinear solutions, including zero-coefficient regression checks. These are numerical tests over specified cases, not experimental validation of a commercial compressor.

The Stage 9B handoff reports 69 enabled regression tests passing at its stated baseline, with one test intentionally disabled. A fresh test count against the latest tree has not been recorded in the reviewed status materials.

### Not yet delivered

| Deliverable | State |
|---|---|
| End-to-end pipe/compressor/pipe example | Not present |
| Transient compressor coupling and coupling convergence study | Not present |
| Reproduction of the hackathon ReLU interpolant in this repo | Not present |
| Qualified surrogate data set, held-out evaluation, and error gate | Not present |
| Surrogate use through the pipeline component interface | Not present |
| CUDA/GPU compressor solver | Not present |
| NASA HECC geometry reconstruction and experimental comparison | Not present |

## Why this has felt open-ended

Most recent work advanced the compressor physics and numerical verification: nonlinear operating envelopes, choking-related studies, frictional diffuser equations, and variable-width wall geometry. That work made the reference model more capable, but it did not itself produce the pipeline coupling or surrogate deliverable.

The NASA HECC data archive was not available in the last recorded handoff (HTTP 403 on 2026-10-13). The repository therefore documents published reference dimensions and numerical checks, but does not claim a reconstructed experimental geometry or data-based validation. This is a model-credibility gap, not a reason to defer a first software-coupling example.

## Immediate work

### Pipeline path

1. Define an adapter from pipe boundary state to the compressor’s inlet state and operating inputs.
2. Return outlet thermodynamic state, power/work, convergence status, and validity information through a matching adapter.
3. Verify one steady pipe–compressor–pipe case with mass and energy checks.
4. Extend that same case to transient coupling and document the coupling iteration or time-step treatment.

### Surrogate path

1. Locate the hackathon’s sample data, code, and exact feature/target conventions.
2. Reproduce its ReLU interpolation result with a fixed, structured holdout.
3. Compare the legacy sample source with the C++ physics model at matched operating conditions.
4. Add physics-model samples with clear provenance and documented validity limits.
5. Qualify a surrogate over a stated domain before exposing it to the pipeline model.

The two paths should share a component interface and test cases. The surrogate can then replace the reference evaluator only where its measured error and domain qualification are acceptable.

## Relevant records

- [Development roadmap](ROADMAP.md)
- [M1.4 nonlinear solver validation](M1.4_VALIDATION.md)
- [M1.5 operating-envelope validation](M1.5_VALIDATION.md)
- [M1.6b.2 incidence validation](M1.6b2_VALIDATION.md)
- [Stage 9B HECC validation handoff](STAGE9B_HECC_HANDOFF.md)
