# TurboCompGPU

TurboCompGPU is a C++20 compressor-modeling project built on
[ThermoGPU](https://github.com/rdm375/ThermoGPU). It is intended to support two related engineering goals:

1. **Pipeline modeling:** provide a compressor component that can be coupled to a transient pipe model.
2. **A semi-realistic surrogate candidate:** continue the ReLU interpolation approach from the original NN Interpolation hackathon, using compressor-operation samples and a physics-based model as possible sources of reference data.

The surrogate is meant to approximate compressor behavior cheaply inside larger studies. The physics solver remains the reference used to generate or check samples; the present model is not yet an experimentally validated compressor station.

## Current status

The default branch contains substantial compressor and diffuser modeling, but the project has not reached either end-to-end goal:

- **Pipeline coupling is not implemented.** The repository has flow-station and compressor operating-point types, but no demonstrated compressor-to-pipeline boundary component or coupled transient case.
- **A trained ReLU surrogate is not implemented here.** The hackathon’s sample-based interpolation is the intended starting point; this repository can contribute physics-based samples, subject to its model limitations.
- **GPU execution is not implemented.** Current solver and validation workflows are CPU based.
- **Experimental compressor validation is pending.** The NASA HECC archive geometry and measurements were not obtained; the archived handoff records an HTTP 403.

See [Project Status](docs/PROJECT_STATUS.md) for the evidence, limitations, and immediate work items, and the [Roadmap](docs/ROADMAP.md) for the two workstreams.

## What is implemented

- Ideal-gas and Peng–Robinson thermodynamics through ThermoGPU.
- Impeller meanline aerodynamics, conservation equations, and nonlinear operating-point solving.
- Compressor-map and operating-envelope calculations with continuation support.
- Impeller and vaneless-diffuser complete-stage calculations, including frictional and variable-width diffuser formulations.
- PCHIP infrastructure for diffuser hub and shroud profiles.
- Incidence-loss experiments with independent algebraic checks.

The incidence study covers 105 paired operating points across three thermodynamic configurations and reports 315 successful nonlinear solutions. Its incidence coefficient is exploratory, not experimentally calibrated. The HECC geometry work provides infrastructure, not a NASA geometry reconstruction or experimental validation.

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

## Project direction

The immediate objective is a small, reviewable pipeline-coupling milestone, developed alongside a reproducible surrogate baseline:

- Define the compressor component inputs and outputs using consistent thermodynamic and flow conventions.
- Demonstrate one compressor between pipe sections in a steady case, then in a transient case.
- Reproduce and qualify the hackathon’s sample-based ReLU interpolation on held-out compressor operation.
- Compare the surrogate against both the available sample source and the physics solver over a clearly defined domain.

See the [Roadmap](docs/ROADMAP.md) for acceptance criteria and sequencing.

## License

GPL-3.0-or-later.
