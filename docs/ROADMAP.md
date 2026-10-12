# TurboCompGPU Development Roadmap

Status: Proposed roadmap, 2026-10-11

## Project direction

Develop a physically validated, high-performance centrifugal
compressor modeling engine suitable for performance-map generation,
accelerator execution, surrogate modeling, and eventual integration
with transient pipeline simulation.

Correctness and physical validity take precedence over acceleration.

## Phase 1 — Physical and numerical foundations

| Milestone | Scope | Status |
|---|---|---|
| M1.0 | Project infrastructure | Complete |
| M1.1 | ThermoGPU 2.0 integration | Complete |
| M1.2 | Thermodynamic foundations | Complete |
| M1.3 | Impeller aerodynamics and conservation | Complete |
| M1.4 | Nonlinear solver, validation, optimization | Complete |
| M1.5 | Operating-envelope exploration | Next |
| M1.6 | Compressor performance maps | Planned |
| M1.7 | Operating limits and validity boundaries | Planned |
| M1.8 | Physical-model extensions and calibration | Planned |

M1.4 numerical and benchmark work is complete.
See M1.4_VALIDATION.md for supporting evidence.

## M1.5 — Operating-envelope exploration

### M1.5a — Operating-point evaluation API

Create a reusable interface accepting a scenario parameter vector
and returning:

- Outlet thermodynamic state.
- Pressure ratio and specific work.
- Shaft power and torque.
- Entropy generation and defined efficiency metrics.
- Nonlinear convergence diagnostics.
- Physical-validity diagnostics.

The interface should be reusable by scalar, parallel, accelerator,
adaptive-sampling, and transient-system workflows.

### M1.5b — Deterministic operating sweeps

Sweep mass flow and rotational speed with fixed inlet conditions
and geometry.

Subsequently introduce inlet pressure and temperature variations.

Preserve complete numerical diagnostics and explicit failure states.

### M1.5c — Continuation methods

Use neighboring converged states to initialize subsequent solves.

Compare continuation against independent initialization.

Measure convergence, iteration counts, robustness, and performance.

Do not assume continuation detects all solution branches.

### M1.5d — Visualization and validation

Generate operating-envelope plots and machine-readable results.

Verify conservation, reference invariance, reproducibility, and
agreement with independently solved operating points.

### M1.5 completion criteria

- Reproducible operating sweeps.
- Explicit distinction between convergence and physical validity.
- Tested operating-point API.
- Documented continuation behavior.
- Saved numerical results and plots.
- Passing regression suite.

## M1.6 — Compressor performance maps

Develop speed lines and performance-map outputs for:

- Pressure ratio.
- Specific work and head.
- Power and torque.
- Efficiency under clearly defined conventions.
- Entropy generation.

Verify map values against the reference solver.

## M1.7 — Operating limits

Investigate:

- Choking and compressibility constraints.
- Surge-related indicators and model limitations.
- Aerodynamic validity boundaries.
- Nonlinear solution multiplicity and branch behavior.

Do not equate nonlinear convergence with physical realizability.
Do not claim predictive surge boundaries without suitable models
and validation.

## M1.8 — Physical-model extensions

Investigate additional loss correlations, diffuser models,
experimental comparison, and parameter calibration.

Separate numerical verification from experimental validation.

## Phase 2 — High-performance execution

| Milestone | Scope |
|---|---|
| M2.0 | Batched operating-point architecture |
| M2.1 | Parallel CPU execution |
| M2.2 | CUDA compressor solver |
| M2.3 | CPU/GPU parity and performance validation |

Measure actual batch throughput before committing to major
accelerator optimizations.

Retain the validated scalar implementation as the numerical oracle.

## Phase 3 — Adaptive compressor-map surrogates

| Milestone | Scope |
|---|---|
| M3.0 | TDAR-guided adaptive sampling |
| M3.1 | Physics-informed coordinate transformations |
| M3.2 | Mixed-basis sparse tensor-product interpolation |
| M3.3 | Error estimation and surrogate validation |

Prioritize compactly supported basis functions.

Candidate one-dimensional families include hierarchical linear
hats and cubic B-splines.

Investigate corrected-flow and corrected-speed coordinates,
as well as normalized operating-envelope coordinates.

Evaluate automatic-differentiation-informed error estimators,
adaptive sample efficiency, and accelerator execution.

Compare surrogate outputs against the validated compressor solver.

## Phase 4 — Transient-system integration

| Milestone | Scope |
|---|---|
| M4.0 | Compressor and transient-pipeline integration |

Investigate coupling compressor performance models to
transient pipeline simulation.

Preserve conservation and consistent thermodynamic conventions.

## Development principles

1. Establish physical and numerical correctness first.
2. Preserve independent validation and regression tests.
3. Distinguish solver failure from physical invalidity.
4. Record provenance for numerical experiments.
5. Benchmark complete workflows, not only isolated kernels.
6. Maintain CPU reference behavior when adding accelerators.
7. Prefer reusable interfaces over backend-specific duplication.

## Roadmap status

M1.5–M4.0 are proposed milestones and may be refined as research
results emerge. Accelerator and surrogate development may proceed
in parallel once suitable reference datasets exist.
