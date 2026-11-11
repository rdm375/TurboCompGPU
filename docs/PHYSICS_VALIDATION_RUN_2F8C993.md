# Physics validation run at `2f8c993`

Run date: 2026-11-11. The repository baseline is remote `main` commit
`2f8c993062b109ccbf652725befea476641c59d4`.

## Decision

The checked-in evidence supports limited air-map screening and numerical
verification of a pure-methane PR thermodynamic path. It does not validate
natural-gas compressor predictions. The full Eckardt A efficiency-map shape
and the independent Choi–Kang mixed-plane efficiency comparison remain poor.

## Fresh results

The Eckardt diffuser-friction executable was compiled directly with `g++`
against the ThermoGPU sources available in the validation environment. CMake
was not installed there, so the full CMake/CTest regression was not run.
Both ideal-gas and PR pseudo-air benchmark runs completed successfully; their
raw output files are in `benchmarks/results-eckardt-impeller-a/remote-2f8c993/`.

| Comparison | Points | Bias | MAE / MAPE | Maximum absolute error |
|---|---:|---:|---:|---:|
| Eckardt A total pressure ratio, ideal gas | 17 | −0.832% | 1.342% MAPE | 3.179% |
| Eckardt A total pressure ratio, PR-air | 17 | −0.802% | 1.328% MAPE | 3.221% |
| Eckardt A total-to-total efficiency, ideal gas | 21 | −3.880 pp | 4.764 pp MAE | 7.661 pp |
| Eckardt A total-to-total efficiency, PR-air | 21 | −3.881 pp | 4.766 pp MAE | 7.663 pp |

For the 14,000 rpm line, measured/model full-line efficiency slopes are
−2.466/+0.943 pp/(kg/s). Excluding the approximate 5.69 kg/s point as a
separate diagnostic gives +3.188/+1.255 pp/(kg/s). The full line remains in
the score; the reduced line is not a replacement acceptance result.

### Independent checks

- **Choi–Kang Fig. 2:** the current qualification script now completes the
  matched-plane comparison. Across three approximate digitized points, PR
  MAPE is 2.49%, mean absolute temperature-rise-coefficient error is 0.0581,
  and mixed-plane efficiency MAE is 6.74 pp. At 3.405 and 4.540 kg/s, model
  efficiencies are low by 8.88 and 8.77 pp. The source digitization tolerance
  is approximately ±1.5 pp. See `choi-kang-fig2/qualification.md`.
- **Eckardt B held-out pressure line:** after fitting an effective slip/head
  parameter on 10,000–14,000 rpm, the held-out 16,000 rpm line has 3.785%
  MAPE and 5.151% maximum absolute error, below the predeclared provisional
  7% maximum-error gate. This is a calibrated, same-family pressure check,
  not a physical slip measurement or natural-gas validation. See
  `eckardt-b-holdout/speedline_holdout.md`.
- **Pure methane thermodynamics:** the tracked 105-point validation passed
  pressure-ratio identities, power identity, baseline comparisons, and
  matched methane PR density identities. This verifies the property/solver
  path under that test; it is not compressor-map validation for a mixture.

### Cross-fluid screens rerun on the same source

The Sandia sCO2 component screen and Schiffmann–Favrat R134a design-point
smoke executable were also compiled directly with `g++` against the same
ThermoGPU source tree. Both executables returned success, but neither result
is a physics-qualification pass:

| Screen | Fresh result | Interpretation |
|---|---|---|
| Sandia sCO2 impeller-only P3 | 12/12 points converged; 4.917% MAPE, 9.684% maximum absolute error | PR/Shomate CO2 properties are not reference quality near the pseudo-critical region; geometry, inlet reconstruction, and slip are provisional; vaned diffuser excluded |
| Sandia CO2 property spot checks | At 305 K, 8 MPa: density −12.33%, Cp +20.56%; at 310 K, 8 MPa: density +1.07%, Cp −12.70% versus the in-source Span–Wagner values | The CO2 property discrepancy alone blocks treating the compressor residual as a loss-model validation |
| Schiffmann–Favrat R134a point | PR model 3.0732 vs measured PR 2.35 (+30.78%); computed efficiency 1.000 vs 0.800 (+20.0 pp) | Smoke executable converges at the ideal efficiency limit but has large PR and efficiency errors; this is not a physical-validation pass |

Raw CSV and console output are under `benchmarks/results-sandia-co2/remote-2f8c993/`
and `benchmarks/results-schiffmann-favrat/remote-2f8c993/`. The Sandia
properties use PR with a NIST Shomate caloric model. The R134a screen uses PR
with constant cp and provisional geometry/slip, so neither should be used to
infer validated real-gas compressor performance.

The NASA HECC handoff remains blocked on archive retrieval. The current NASA
NTRS guide labels the archive public and gives its URL, but a fresh HTTP HEAD
request in this environment still returned 403 AccessDenied. The NASA paper
describes open geometry, performance, and detailed aerodynamic data, but no
HECC measurements were ingested or scored in this run. HECC is an air/turbine
compressor dataset, so even successful use would strengthen aerodynamic
validation rather than close the natural-gas gate.

## Why PR-air barely changes this benchmark

The source selects either ideal gas or a Peng–Robinson pseudo-pure dry-air
model. The PR-air model itself is configured with the same constant
`cp=1004.5 J/(kg K)` as ideal air; the comparison therefore isolates the EOS
pressure-volume nonideality and associated state inversions rather than
variable heat capacity or a natural-gas mixture.

Compressibility is present in the state calculations: the EOS provides density
and thermodynamic states, and the selected model is used for entropy/enthalpy
inversions and isentropic sound speed in the choking calculation. The loss
correlations do not generally include an explicit compressibility-factor
correction. They are predominantly velocity/geometry correlations, with
Reynolds number and density ratios entering selected terms. The EOS therefore
affects them indirectly through density, state, and choking inputs.

At the low absolute pressures and moderate temperatures in this air map,
PR-air is close to ideal. The maximum absolute change in the smooth-wall
predictions between the fresh ideal-gas and PR-air runs is 0.000874 in stage
pressure ratio and 0.000112 in efficiency. The corresponding validation errors
are essentially unchanged. This is expected for this air condition and does
not imply that real-gas effects are negligible for high-pressure natural gas.
The current target-specific natural-gas conditions and composition have not
been run through this Eckardt benchmark.

## Reproducibility issue found

The tracked PR-air `loss_term_ablation.csv` predates the current C++ source:
it contains 51 scenarios, including 20 dense `flow_mix_sweep_chi_*` cases,
while the current source emits 31 scenarios. The fresh run outputs here have
31 scenarios. Refresh or remove the stale tracked artifacts before treating
all stored ablation output as reproducible evidence.

The Choi–Kang qualification script previously aborted if the dense mixing
sweep was absent. Its updated behavior is to write the matched-plane results
and explicitly report that the optional mixing-flow diagnostic was not run.
This change does not alter the physical model or benchmark calculations.

## Natural-gas paper comparison

Two open-access papers and the five Chen et al. (2025) operating points are
now archived and compared with the current evidence in
[`NATURAL_GAS_PAPER_COMPARISON.md`](NATURAL_GAS_PAPER_COMPARISON.md). This is a
source-data and methods comparison only: current compressor runs are on the
Eckardt air machine, so there are no like-for-like TurboCompGPU residuals on
either natural-gas machine.

## Scope needed to close natural-gas validation

A blind, geometry-matched natural-gas compressor comparison is still needed.
It must include measured gas composition, machine geometry, station
conventions, multiple speeds/flows, relevant performance measurements and
uncertainties, plus predeclared output-specific acceptance limits. Until then,
report natural-gas predictive accuracy as unqualified.
