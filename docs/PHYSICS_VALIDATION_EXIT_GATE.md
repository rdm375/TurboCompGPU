# Physics validation exit gate

Reviewed: 2026-11-11

## Decision

**Physics validation for natural-gas compressor predictions is not complete.**
The current Eckardt work is useful aerodynamic evidence for air. It does not
qualify the model for a natural-gas compressor, and the latest benchmark
selects only `ideal-gas` or `pr-air`.

The broader repository does have a passing pure-methane Peng–Robinson
thermodynamic comparison (105 operating points). That qualifies the methane
property and solver path covered by that benchmark; it does not qualify a
multicomponent pipeline-gas mixture or validate the meanline loss model on a
natural-gas compressor map.

Accelerator implementation can be evaluated against a frozen CPU reference
once the scalar reference is reproducible. Such results would establish
backend agreement and speed, not natural-gas predictive accuracy. If the
project requires physical validation to precede any accelerator work, the
natural-gas gate below remains a blocker.

## What the existing evidence establishes

| Evidence | Conclusion | Limit |
|---|---|---|
| Solver, conservation, derivative, and property-path tests | Numerical implementation has meaningful verification coverage | A fresh full regression on the current tree is not recorded here; numerical verification is not experiment agreement |
| Eckardt A and Choi–Kang comparisons | The meanline stage has been screened against same-family air data and a few station-matched plotted points | Some inputs are digitized; Choi–Kang is one speed line; neither dataset is natural gas |
| Published Set 2042 scenario | A published correlation bundle can reproduce pressure-ratio trends reasonably at selected Eckardt A points | The implementation uses proxies and does not reproduce the source paper's full geometry or algorithm; efficiency remains discrepant |
| Loss, mixing, energy-partition, and slope sweeps | The residual is sensitive to loss closure and flow; energy bookkeeping can be inspected | A closer fit to the development map is not independent validation or proof that a particular correlation is correct |
| Real-gas Eckardt run | Peng–Robinson air gives a useful check of the EOS path at air conditions | This is not a natural-gas mixture test; `pr-air` is explicitly the air composition |

## Fresh checks on 2026-11-11

The stored thermodynamic comparison was rerun and passed: 105 operating
points, the M1.6a baseline comparison, pressure-ratio identities, the power
identity, and matched methane EOS density identities. The reported mean
Peng–Robinson-minus-ideal methane differences were 0.00437 in total pressure
ratio and 130.8 J/kg in specific work. These are comparisons between the two
models, not prediction errors against a compressor test map.

The current Eckardt source was compiled directly with `g++` against the local
ThermoGPU sources because CMake is unavailable in this checkout. The
Peng–Robinson-air map produced 17 paired pressure-ratio rows and 21 efficiency
rows for the smooth-wall, Galvas-incidence case:

| Quantity | Points | Bias | MAE | Maximum absolute error |
|---|---:|---:|---:|---:|
| Total pressure ratio | 17 | −0.802% | 1.328% | 3.221% |
| Total-to-total efficiency | 21 | −3.881 pp | 4.766 pp | 7.663 pp |

At 14,000 rpm the measured efficiency slope across all five digitized points
is −2.466 pp/(kg/s), while the current model gives +0.942 pp/(kg/s). On the
four ordinary-flow points below 5.69 kg/s, the slopes are +3.188 and
+1.255 pp/(kg/s), respectively. The high-flow point remains part of the
full-map score; the four-point result is a separate diagnostic. This means
pressure-ratio agreement is useful, but the efficiency map shape is not
validated.

The first run exposed a zero-loss limit bug: the entropy inversion produced a
small negative roundoff value when `internal_loss == 0`, causing the
`without_impeller_losses` diagnostic at 10,000 rpm and 2.34 kg/s to fail a
nonnegative guard. The closure now assigns exactly zero entropy generation for
that exact zero-loss limit. The rerun completed all ablations and wrote the
full diagnostic CSV. This is a numerical edge-case correction, not a change to
the nonzero-loss physics. The 14,000 rpm efficiency-slope mismatch remains.
The complete run artifacts and pointwise slope tables are under
`benchmarks/results-eckardt-impeller-a/physics-validation-current/pr-air/`.

The approximate 5.69 kg/s Eckardt A efficiency point is a separate high-flow
stress point. It must not be allowed to dominate the ordinary-flow assessment
or be silently discarded from the final report.

## Minimum natural-gas gate

Before claiming the model is physically validated for a target natural-gas
compressor, freeze one test article and obtain:

1. **Gas:** a measured mole composition (including treatment of C5+), the
   source EOS/property method, and operating-condition ranges.
2. **Machine:** impeller and diffuser geometry sufficient to instantiate this
   meanline model, with station and total/static conventions.
3. **Measurements:** multiple speeds and flows with mass or volume flow,
   inlet/outlet total pressure and temperature, shaft power or torque where
   available, and source uncertainty.
4. **Blind comparison:** freeze equations, loss settings, and geometry mapping
   before scoring at least one held-out speed line or operating range.
5. **Acceptance:** declare output-specific residual limits from combined
   measurement, input, and digitization uncertainty before revealing the
   holdout. Report bias, MAE/RMSE, maximum error, trend with flow/speed, and
   solver failures; do not collapse the outputs into one score.

Without machine geometry, field data can check a station-level or map
conversion model, but it cannot validate this geometry-based aerodynamic
solver. Without a specified gas composition, Peng–Robinson air is not a
substitute for natural gas.

## Practical accelerator entry gate

The accelerator track needs a **frozen scalar CPU oracle**, not a claim that
all compressor physics is settled. Before using accelerator results for
optimization, require:

- a recorded source revision, compiler/build settings, gas composition,
  geometry, and fixed operating-point set;
- current CPU regression and conservation results;
- CPU outputs for pressure ratio, efficiency, work, and solver status at
  interior and boundary points;
- accelerator outputs compared pointwise against that CPU oracle with
  predeclared tolerances, including derivatives if optimization uses them;
- no accelerator-only changes to correlations or physical constants.

This allows backend parity and throughput experiments to be judged cleanly.
It does not lift the natural-gas validation limitation above.

## Natural-gas dataset lead and limit

An open 2025 *Energies* paper reports five actual gas-storage compressor
operating points and natural-gas composition, alongside pressure ratio and
power comparisons. The authors validate a performance-curve conversion
method; the public paper does not provide the matched detailed geometry and
uncertainty needed to score this meanline solver as a blind aerodynamic
prediction. Treat it as a possible thermodynamic/map-level regression lead,
not as a completed qualification dataset.

Source: H. Chen et al., “The Influence of Variable Operating Conditions and
Components on the Performance of Centrifugal Compressors in Natural Gas
Storage Reservoirs,” *Energies* 18(15), 3930 (2025),
[doi:10.3390/en18153930](https://doi.org/10.3390/en18153930). The source
describes its method as similarity-based performance-curve conversion and
states that its heat-capacity ratio is held constant; its five field points
are not a geometry-resolved validation of this code.

A more directly relevant paper is Y. Lin et al., “Healthy-State Performance
Modeling of a Multistage Natural Gas Centrifugal Compressor Using a
CFD-Generated Baseline and Factory-Data Correction,” *Processes* 14(13), 2154
(2026), [doi:10.3390/pr14132154](https://doi.org/10.3390/pr14132154). It
describes a three-stage pipeline compressor, a methane-dominant mixture, PR
checked against GERG-2008, and comparisons against factory data. The authors
report that their raw CFD baseline has up to 9.802% polytropic-efficiency
error on an independent 44-point subset, and that a factory-data residual
correction reduces the maximum to 0.899%. The article explicitly says the
data are confidential; its reconstructed geometry and pointwise measurements
are not available for an independent run here. It is evidence that even a
validated CFD baseline can retain systematic map residuals, not an available
qualification dataset for this meanline model.

## Closeout rule

Mark **air meanline validation** complete only after the ablation runner exits
cleanly, the current source tree passes its recorded regression, and the
independent air-map report is frozen with uncertainty limits and an explicit
disposition of the efficiency-slope residual. Mark **natural-gas compressor validation**
complete only after the blind, geometry-matched mixture dataset above passes
its predeclared acceptance criteria. Until then, report the model as a
numerically verified research reference with limited air-map screening and
unqualified natural-gas predictive accuracy.
