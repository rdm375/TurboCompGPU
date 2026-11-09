# TurboCompGPU

TurboCompGPU is a C++20 compressor-modeling project built on
[ThermoGPU](https://github.com/rdm375/ThermoGPU). It is intended to support two related engineering goals:

1. **Pipeline modeling:** provide a compressor component that can be coupled to a transient pipe model.
2. **A semi-realistic surrogate candidate:** continue the ReLU interpolation approach from the original NN Interpolation hackathon, using compressor-operation samples and a physics-based model as possible sources of reference data.

The surrogate is meant to approximate compressor behavior cheaply inside larger studies. The physics solver remains the reference used to generate or check samples; the present model is not yet an experimentally validated compressor station.

## Current status

The default branch contains substantial compressor and diffuser modeling. The initial software pipe-coupling and surrogate workflows are present, but neither is fully qualified for pipeline use:

- **Pipe/compressor coupling is implemented and is isothermal/quasi-steady.** The optional `pipe_sim` adapter connects upstream and downstream pipe segments, supports recycle accounting, and includes single-step and multi-step checks. A Dell run of the new recycle disturbance case and time-step refinement are still pending.
- **The ReLU surrogate training/export and held-out-speed-line evaluator are included.** They approximate current C++-generated samples, not the unavailable hackathon dataset or experimental truth. No current holdout report or pipeline-use error gate is recorded.
- **Guide-vane mapping has a calibrated callback interface.** A selected-machine physical map is still missing. A public TUM experiment is a cross-machine digitization candidate, not calibration for this compressor.
- **Surge protection requires an external map.** The current solver does not predict surge; keep anti-surge unavailable until a selected-machine map is independently qualified.
- **GPU execution is not implemented.** Current solver and validation workflows are CPU based.
- **Experimental compressor validation is partial.** Eckardt, Choi–Kang, and Sandia screens are documented with limitations. The NASA HECC archive geometry and measurements were not obtained; the archived handoff records an HTTP 403.

See [Project Status](docs/PROJECT_STATUS.md) for the evidence, limitations, and immediate work items, and the [Roadmap](docs/ROADMAP.md) for the two workstreams.
The current held-out-surrogate, recycle-transient, and guide-vane/surge data
record is in [Project gap-close data and run record](docs/PROJECT_GAP_CLOSE_2026-11.md).

## What is implemented

- Ideal-gas and Peng–Robinson thermodynamics through ThermoGPU.
- Impeller meanline aerodynamics, conservation equations, and nonlinear operating-point solving.
- Compressor-map and operating-envelope calculations with continuation support.
- Impeller and vaneless-diffuser complete-stage calculations, including frictional and variable-width diffuser formulations.
- PCHIP infrastructure for diffuser hub and shroud profiles.
- Incidence-loss experiments with independent algebraic checks.

The incidence study covers 105 paired operating points across three thermodynamic configurations and reports 315 successful nonlinear solutions. Its incidence coefficient is exploratory, not experimentally calibrated. The HECC geometry work provides infrastructure, not a NASA geometry reconstruction or experimental validation. For exact current limits and next work, see [Project Status](docs/PROJECT_STATUS.md) and the [Roadmap](docs/ROADMAP.md).

## Validation and performance

The repository contains extensive numerical regression, manufactured-solution, convergence, sensitivity, and independent-residual checks. These establish numerical behavior for the tested cases; they do not establish a predictive surge boundary or experimental accuracy.

The optimized CPU residual implementation was measured at roughly 48–61× lower median execution time than the original full-residual implementation on an Intel Core i7-6920HQ. Including independent validation, the reported reduction was about 12–16×. These are historical CPU comparisons, not GPU speedups or compressor-model accuracy claims.

Milestone reports:

- [M1.4 nonlinear solver validation](docs/M1.4_VALIDATION.md)
- [M1.5 operating-envelope validation](docs/M1.5_VALIDATION.md)
- [M1.6b.2 incidence validation](docs/M1.6b2_VALIDATION.md)
- [Stage 9B HECC validation handoff](docs/STAGE9B_HECC_HANDOFF.md)

## Build

Requirements: C++20 compiler, CMake 3.20+, and an installed ThermoGPU 2.0 package.

    cmake -S . -B build-cpu \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=/path/to/ThermoGPU/install
    cmake --build build-cpu -j2
    ctest --test-dir build-cpu --output-on-failure

The exact test set depends on the configured build options. Consult CMake output and the milestone reports when comparing test counts.

## Benchmark

    ./build-cpu/turbocompgpu_nonlinear_wall_clock \
      benchmarks/results-wall-clock/timings.csv 1000 100
    python benchmarks/analyze_wall_clock.py \
      benchmarks/results-wall-clock/timings.csv

## ReLU surrogate baseline

A first training script uses the repository’s existing compressor-map CSV to fit separate ReLU networks for IdealGas and PengRobinson. It predicts total outlet temperature, total outlet pressure, and specific work from mass flow and speed under the benchmark’s fixed inlet and geometry conditions. The data contain 35 unique operating conditions per thermodynamic model; independent and continuation rows are deduplicated.

Run it with:

    python3 benchmarks/train_relu_surrogate.py \
      benchmarks/results-compressor-map/map.csv \
      --output-dir build/relu-surrogate

This is a model-generated interpolation baseline, not an experimental or legacy-station surrogate. See [ReLU surrogate baseline](docs/RELU_SURROGATE_BASELINE.md) for the domain and structured speed-line holdout method.

## Project direction

The remaining project work is to run and review the surrogate holdout, verify the multi-step recycle transient on Dell, obtain selected-machine vane/recycle/surge data, and validate against measured station transients:

See the [Roadmap](docs/ROADMAP.md) for acceptance criteria and sequencing.

## License

GPL-3.0-or-later.
