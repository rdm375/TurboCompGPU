# TurboCompGPU

Physics-based centrifugal compressor modeling with real-gas
thermodynamics and high-performance numerical methods.

TurboCompGPU is an open-source scientific computing project
for modeling centrifugal compressor performance using
ThermoGPU 2.0, meanline aerodynamics, and validated nonlinear
numerical methods.

The project is developed in C++20 and licensed under
GPL-3.0-or-later.

## M1.6b.2 — Incidence-aware compressor modeling

The current development milestone introduces a physics-based,
flow-dependent incidence-loss closure into the nonlinear
centrifugal compressor solver.

Validation includes:

- Three thermodynamic models: air-like ideal gas,
  methane ideal gas, and methane Peng-Robinson.
- 105 paired compressor operating points.
- 315 successful nonlinear solutions, including
  zero-coefficient regression checks.
- Analytical validation of the zero-incidence locus.
- Independent algebraic incidence-loss validation.
- Static and total pressure-ratio comparisons.

The incidence coefficient is exploratory, not experimentally
calibrated. The current cross-gas comparison uses separately
calibrated blade angles.

GPU execution remains planned rather than implemented.

See [M1.6b.2 validation](docs/M1.6b2_VALIDATION.md)
for the governing equations, results, limitations, and
reproducibility instructions.

## Current status

**M1.4 complete: Validated nonlinear compressor solver
and optimized CPU implementation.**

The current implementation supports:

- Ideal-gas and Peng-Robinson thermodynamics.
- Thermodynamic reference-invariant calculations.
- Impeller meanline geometry and aerodynamics.
- Mass, angular-momentum, and energy conservation.
- Slip and entropy-generation closures.
- Nonlinear operating-point solution.
- Optimized thermodynamic residual evaluation.
- Independent physical validation.
- Deterministic failure diagnostics.
- Reproducible CPU benchmarking.

GPU execution is planned but has not yet been implemented.

## Performance

The optimized nonlinear solver achieves approximately
**48–61x lower median execution time** than the original
full-residual implementation across six benchmark operating
conditions.

Including independent physical validation, acceleration
remains approximately **12–16x**.

Measurements were collected on an Intel Core i7-6920HQ
using GCC 15.2.0 with -O3 optimization.

| Thermodynamic model | Fast solve | Fast + validation |
|---|---:|---:|
| Ideal gas | 47.9–49.7x | 11.9–12.5x |
| Peng-Robinson | 61.0–61.1x | 15.9x |

These are CPU wall-clock measurements, not GPU speedups.

See [M1.4 Validation](docs/M1.4_VALIDATION.md)
for numerical results, methodology, limitations, and
reproducibility information.

## Numerical validation

The M1.4 regression suite contains 26 CTest tests.

Additional experiments include:

- 54 thermodynamic reference-shifted solves.
- 384 paired convergence-basin cases.
- Deterministic coverage of all seven nonlinear
  termination statuses.
- Independent five-residual physical validation.

The full and optimized solvers exhibited no convergence
diagnostic disagreements in the tested convergence basin.

## Dependencies

- C++20 compiler.
- CMake 3.20 or newer.
- ThermoGPU 2.0.

Python dependencies for benchmark analysis are listed in
requirements.txt.

## Build

    cmake -S . -B build-cpu \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH=/path/to/ThermoGPU/install

    cmake --build build-cpu -j2

## Test

    ctest --test-dir build-cpu --output-on-failure

## Benchmark

    ./build-cpu/turbocompgpu_nonlinear_wall_clock \
        benchmarks/results-wall-clock/timings.csv \
        1000 100

    python benchmarks/analyze_wall_clock.py \
        benchmarks/results-wall-clock/timings.csv

## Roadmap

### Completed

- Thermodynamic integration.
- Impeller meanline modeling.
- Conservation and entropy-loss closures.
- Nonlinear operating-point solving.
- Optimized CPU residual evaluation.
- Numerical validation and benchmarking.

### Next: M1.5

- Compressor operating-envelope exploration.
- Mass-flow and rotational-speed sweeps.
- Compressor performance-map generation.
- Efficiency and pressure-ratio characterization.
- Physical operating-limit investigation.

### Future research

- OpenMP and CUDA acceleration.
- Batched GPU operating-point evaluation.
- Geometry- and physics-informed compressor-map surrogates.
- Adaptive sampling and sparse tensor-product interpolation.
- Integration with transient pipeline simulation.

See [Development Roadmap](docs/ROADMAP.md) for milestone details.

## Development philosophy

1. Establish a correct scalar numerical reference.
2. Validate against analytical and independent references.
3. Implement and differentially validate parallel backends.
4. Characterize performance using reproducible benchmarks.

## ThermoGPU

TurboCompGPU builds upon:

https://github.com/rdm375/ThermoGPU

## License

GPL-3.0-or-later.


## M1.5 — Operating-envelope validation

The M1.5 operating-envelope implementation evaluates centrifugal
compressor operating points over deterministic mass-flow and
rotational-speed grids.

The validated experiment includes:

- 70 operating conditions across IdealGas and PengRobinson models.
- 140 nonlinear solves using independent and continuation initialization.
- 140/140 converged solutions.
- 140/140 independently validated solutions.
- 140/140 solutions passing basic physical-consistency checks.
- Continuation improvements at all 60 points with available history.
- Eight operating-envelope visualizations.

Continuation reduced mean Newton iterations by 34.78% for IdealGas
and 42.45% for PengRobinson.

These are numerical iteration reductions, not measured wall-clock
speedups.

Surge and choke limits remain unevaluated. The present model uses
prescribed entropy generation rather than a predictive aerodynamic
loss model.

See [M1.5 validation](docs/M1.5_VALIDATION.md) for the complete
experiment and reproducibility instructions.

