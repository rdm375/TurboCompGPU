# Physics validation exit gate

Reviewed: 2026-11-11

## Decision

**Physics validation for natural-gas compressor predictions is not complete.**
The current Eckardt work is useful aerodynamic evidence for air. It does not
qualify the model for a natural-gas compressor. The benchmark now also exposes
the Lin et al. fixed-composition PR natural-gas screening mixture.

The broader repository also has a passing pure-methane Peng–Robinson
thermodynamic comparison (105 operating points). The new mixture property-grid
run is a separate code-path check and does not qualify the meanline loss model
on a natural-gas compressor map.

ThermoGPU now has a fixed-composition `thermomodel::PengRobinsonMixture` that
implements the stage-facing `Model` and `DifferentiableModel` interfaces,
including mixture density, enthalpy, entropy, sound speed, and analytic
T,p-derivatives. A reproducible eight-species composition from Lin et al.
(2026), Table 1, is wired into the TurboCompGPU property-grid benchmark and
available as the `pr-natural-gas-lin2026` stage benchmark selector. This closes
the prior software interface gap. The selected benchmark composition uses zero
binary interaction parameters and constant component ideal-gas heat capacities
because the paper does not publish those inputs; it remains a screening setup,
not an independently reference-qualified natural-gas property package.

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
| Lin et al. natural-gas PR property grid and NIST GERG-2008 comparison | The fixed-composition PR model evaluates the 20 published inlet T,p combinations and has a reproducible independent reference comparison | Maximum errors are 2.914% in Z, 3.004% in density, and 2.251% in cp/cv; these exceed Lin's reported maxima. The 19.512% common-reference enthalpy percentage is not directly comparable because Lin's enthalpy datum is unspecified. The reference-invariant isothermal pressure-span enthalpy-change error is 5.9–7.1%. Zero kij and constant component cp remain screening assumptions. Gate 1 fails. |
| Natural-gas stage interface smoke | `pr-natural-gas-lin2026` compiled and completed synthetic 14,000 rpm stage points at 0.25–2.00 kg/s | Synthetic points and Eckardt geometry are only an interface smoke; applying the air-map flows to this gas caused existing physical checks to reject a point, as expected for an unmatched map |
| Sandia CO2 and Schiffmann–Favrat cross-fluid screens | Additional EOS/component paths compile and converge | Sandia PR/Shomate spot properties differ from Span–Wagner by up to 20.6% in Cp; Schiffmann–Favrat returns eta=1.000 and +30.8% PR error; neither is a qualification pass |

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

The 2026-11-11 Sandia sCO2 rerun converged at 12/12 points with 4.917% P3
MAPE and 9.684% maximum error, but uses PR plus a Shomate caloric model and
provisional geometry/blockage/slip. Its 305 K, 8 MPa reference spot check has
−12.33% density error and +20.56% Cp error. The 2026-11-11 Schiffmann–Favrat
R134a smoke run reports PR 3.0732 against 2.35 (+30.78%) and efficiency 1.000
against 0.800 (+20.0 pp). The efficiency is at the ideal limit, while the
pressure and efficiency errors are large. Both executables exit successfully
on convergence; the reported values make clear that convergence is not
physical acceptance.
Artifacts are under the `remote-2f8c993` subdirectories in
`benchmarks/results-sandia-co2/` and `benchmarks/results-schiffmann-favrat/`.

The official NASA HECC guide currently identifies its archive as public and
provides the archive URL, and its related NASA presentation describes public
geometry, performance, and detailed aerodynamic data. However, a fresh request
to the archive URL from this validation environment returned HTTP 403
AccessDenied, so the archive was not ingested or evaluated. HECC would provide
an independent air-compressor benchmark, not a natural-gas validation set.

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
4. **Thermodynamic mixture path:** verify the implemented fixed-composition
   PR model against GERG-2008 or another reference-quality backend over the
   target compressor's pressure/temperature envelope. Replace screening cp
   constants with documented caloric correlations and use published or
   independently justified binary interaction parameters.
5. **Blind comparison:** freeze equations, loss settings, and geometry mapping
   before scoring at least one held-out speed line or operating range.
6. **Acceptance:** declare output-specific residual limits from combined
   measurement, input, and digitization uncertainty before revealing the
   holdout. Report bias, MAE/RMSE, maximum error, trend with flow/speed, and
   solver failures; do not collapse the outputs into one score.

Without machine geometry, field data can check a station-level or map
conversion model, but it cannot validate this geometry-based aerodynamic
solver. The Lin mixture property run does not validate Chen et al.'s different gas
composition or the compressor aerodynamics. Peng–Robinson air remains no
substitute for natural gas.

### Current three-gate status (2026-11-11)

| Gate | Status | Evidence / closure condition |
|---|---|---|
| Mixture EOS versus GERG/REFPROP | **Fail** | The NIST GERG-2008 screen exceeds the reported Lin maxima for Z, density, and cp/cv. The enthalpy percent cannot be compared directly because the paper's datum is unspecified; reference-invariant enthalpy-change diagnostics also show a 5.9–7.1% discrepancy. See `benchmarks/results-natural-gas/lin2026_pr_vs_gerg2008_summary.md`. Determine whether sourced mixture/caloric parameters resolve the mismatch, or use a validated reference-quality mixture EOS. Do not tune to aggregate maxima alone. |
| Compressor prediction versus matched natural-gas measurements | **Blocked** | Chen et al. publish field/map-conversion points but not the stage geometry and flow reference basis needed by this solver; Lin et al. withhold pointwise machine data and geometry. A matched machine dataset with composition, geometry, station definitions, flow basis, inlet/outlet states, power, and uncertainties is required. |
| CPU versus accelerator numerical agreement | **Not run here** | This checkout has no CMake, CUDA compiler, NVIDIA device nodes, or accelerator runtime. Run a frozen CPU oracle and pointwise CUDA/accelerator comparisons on the project machines; compare statuses and outputs (and derivatives if optimization uses them) over interior and boundary cases. |

The repository does not currently have a credible all-three-gates pass. These
states are explicit failures/blockers, not waived gates.

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

## Natural-gas paper comparison

The repo now archives two open-access natural-gas compressor papers and a
source-labeled transcription of the five Chen et al. (2025) gas-storage points.
The row-by-row published values, gas composition, and comparison against the
current TurboCompGPU evidence are documented in
[`NATURAL_GAS_PAPER_COMPARISON.md`](NATURAL_GAS_PAPER_COMPARISON.md). The Chen
points are useful for a later map-level mixture regression, but they do not
provide matched geometry or uncertainty for blind meanline validation. Lin et
al. (2026) report aggregate PR/GERG and compressor errors, while explicitly
withholding pointwise data and geometry. No direct TurboCompGPU-versus-paper
natural-gas residual has been computed.

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

## Additional non-NASA public-data survey (2026-11-11)

The Dewar et al. (2019) open-access paper, tabulated 9-point map, and 12 tap
radii are in the repository. The diffuser pressure measurements are plotted in
Figure 8 and still need calibrated digitization. Rossbach et al. (2015) is
archived under its CC BY-NC-ND 4.0 license; approximate experimental Figure 9
pressure-recovery points are transcribed and their qualitative trends pass a
small analysis script. The Purdue Methel thesis is catalogued but its download
was blocked here, so no values are transcribed. DATED remains a 1D model-output
dataset, not independent experimental evidence.

These references broaden the available aerodynamic evidence, but do not yet
close an external model-comparison gate. The current Eckardt map executable
uses Eckardt A geometry and inlet assumptions, so scoring Dewar, Purdue, or
RWTH points through it would be an invalid cross-machine comparison. The
Rossbach figure read-off is only a source-data trend check. A case-driven
geometry/station interface and source-specific measurement operators are
required before model scores are valid. None of these air cases substitutes
for a matched natural-gas stage dataset.
