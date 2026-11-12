# Compressor method and model-reference audit

Reviewed: 2026-11-11  
Code baseline: `cfc5ad8d32dfe9a74f0abbf00723e6f1ea0bf37b`

## Purpose and status

This audit connects the equations and reported outputs in the current code to
primary method and model references. A reference makes an equation
traceable; it does not establish that the code implements the full reference
method or that the compressor prediction is accurate.

The new reference set closes a **bibliography and method-definition gap**. It
does not close the remaining validation gates. The current Eckardt benchmark
reports total-to-total efficiency, not PTC-10 polytropic efficiency. No
PTC-10 numerical-method comparison, full measurement-uncertainty budget, or
matched natural-gas stage-map validation has been completed.

## Equation-to-code map

| Quantity or model | Reference basis | Current implementation | Audit result / limit |
|---|---|---|---|
| Impeller Euler work | Euler turbomachinery equation; Aungier and Dixon–Hall | The core impeller uses its velocity-triangle/Euler relation in [`impeller_system.hpp`](../include/turbocompgpu/impeller_system.hpp); the Eckardt runner records Euler work and stage enthalpy rise in [`eckardt_diffuser_friction_sweep.cpp`](../benchmarks/eckardt_diffuser_friction_sweep.cpp). | The internal residuals check bookkeeping and conservation. They do not prove that the reconstructed geometry, slip, or loss partition represents Eckardt A. |
| Slip factor | Wiesner (1967), with Busemann comparison; Aungier and Dixon–Hall | The Eckardt map initializer uses `1 - sqrt(cos(30 deg))/Z^0.7` with 20 blades in [`eckardt_diffuser_friction_sweep.cpp`](../benchmarks/eckardt_diffuser_friction_sweep.cpp). | **Convention audit still required.** Verify that the 30° input is the angle convention required by the cited expression, including any complementary-angle conversion and applicable radius-ratio correction. The code formula is documented; its equivalence to the selected published form has not yet been numerically qualified against an independent slip reference. |
| Station thermodynamics and total states | ASME PTC 10-2022; ISO 5389:2005; property-model documentation | Stagnation states are formed from local flow stations with the selected property model. The Eckardt benchmark exposes `ideal-gas`, `pr-air`, and `pr-natural-gas-lin2026` paths. | Natural-gas PR mixture checks against GERG currently fail published Lin bounds; a selected EOS is not itself a property validation pass. |
| Total-to-total efficiency | ASME PTC 10-2022 test-result definitions; compressor thermodynamics | The Eckardt benchmark computes `eta_tt = (h03s-h01)/(h03-h01)`, where `s03s=s01` at predicted `p03`. This is explicit in the benchmark implementation. | This is an isentropic total-to-total efficiency. Do not compare it directly with polytropic efficiency or a paper's metric until station and efficiency definitions match. |
| Polytropic head and efficiency | ASME PTC 10-2022; Taher–Evans 2021; Srinivasan et al. 2024 | The Eckardt runner now integrates a constant-ηp EOS path in pressure coordinates using analytic (dh/dp|_T) and (dh/dT|_p), then shoots on endpoint temperature. It reports `predicted_eta_p` and integrated polytropic head separately from `predicted_eta_tt`. | Ideal-gas and pressure-dependent nonideal-gas analytic checks pass. This is an independently derived diagnostic, **not a claimed PTC 10 or Taher–Evans implementation**. A published numerical example has not yet been reproduced with matched source conventions and EOS. |
| Impeller loss correlations | Aungier (2000); Kovář et al. (2021); Yang et al. (2023) | Correlation bundles are implemented in the Eckardt diagnostic runner and described in [`ECKARDT_PUBLISHED_LOSS_SETS.md`](ECKARDT_PUBLISHED_LOSS_SETS.md). | Selected Set 2042/2137 formulas do not reproduce the source's full geometry or Galvas algorithm. Documented throat and diffuser proxies remain. These are screening alternatives, not a validated closure. |
| Vaneless diffuser | Japikse (1996); Dixon–Hall (2014); Avramenko et al. (2015) | A variable-width diffuser integrates mass, angular momentum, energy, and entropy equations; Avramenko observations are transcribed in [`avramenko_etal_2015_diffuser.csv`](../data/validation/avramenko_etal_2015_diffuser.csv). | The Avramenko rows are section-average **static** pressure ratios and are not yet scored by a geometry-matched diffuser run. Do not compare them to total-pressure ratio. |
| Test uncertainty and decision limits | ASME PTC 10-2022; ASME PTC 19.1-2018 (R2024); ISO 5389:2005 | Reports give point residuals and map-level bias/MAE; they do not propagate sensor, gas-property, geometry, and digitization uncertainty into output confidence intervals. | **Open method gap.** A pass/fail threshold should be based on combined uncertainty and declared before evaluating a holdout; current digitized Eckardt scores are descriptive. |

## Method references and intended use

1. **ASME PTC 10-2022, _Axial and Centrifugal Compressors_.** Use as the
   controlling source for compressor test preparation, instrumentation,
   thermodynamic property handling, correction to specified conditions,
   result definitions, and uncertainty rules. It covers test procedures and
   output quantities; it is not a loss-correlation handbook. [Official ASME
   record](https://www.asme.org/codes-standards/find-codes-standards/axial-and-centrifugal-compressors).
2. **ISO 5389:2005, _Turbocompressors — Performance test code_.** Independent
   test-method reference with broad turbocompressor scope. ISO records the
   edition as confirmed in 2022 and current. [Official ISO
   record](https://www.iso.org/standard/39665.html).
3. **ASME PTC 19.1-2018 (R2024), _Test Uncertainty_.** Use to build and
   propagate Type A/Type B and random/systematic contributions. [Official
   ASME record](https://www.asme.org/codes-standards/find-codes-standards/test-uncertainty).
4. **Taher and Evans (2021), _Centrifugal Compressor Polytropic
   Performance—Improved Rapid Calculation Results—Cubic Polynomial Methods_.**
   Provides a real-gas temperature–entropy path construction and a published
   multisegment calculation. It is a method candidate and numerical
   cross-check, not the text of PTC 10. [DOI
   10.3390/ijtpp6020015](https://doi.org/10.3390/ijtpp6020015).
5. **Aungier (2000), _Centrifugal Compressors: A Strategy for Aerodynamic
   Design and Analysis_; Japikse (1996), _Centrifugal Compressor Design and
   Performance_; Dixon and Hall (2014), _Fluid Mechanics and Thermodynamics of
   Turbomachinery_.** Use as complementary design/model references for meanline
   assumptions, velocity triangles, impeller/diffuser behavior, and correlation
   context. These books are not experimental scorecards for this code.

## Efficiency definitions that must stay separate

For a real gas, the Gibbs relation along a reversible reference path is

\[
dh = T\,ds + v\,dp.
\]

For a compression path with positive enthalpy rise, polytropic efficiency can
be expressed as the ideal differential work divided by actual enthalpy rise,
\(\eta_p=\int v\,dp/\Delta h\), equivalently
\(\eta_p=1-\int T\,ds/\Delta h\), when the path and sign convention are
consistent. Taher–Evans approximate the real-gas \(T(s)\) path with linked
cubic segments. For segment \(i\), their Eq. 44 approximates the degraded
energy integral as

\[
I_i=\int_{s_i}^{s_{i+1}}T(s)\,ds \approx
\frac{T_i+T_{i+1}}{2}\,\Delta s_i
-\frac{E_{i+1}-E_i}{12}(\Delta s_i)^2,
\qquad E_i=\left.\frac{dT}{ds}\right|_{s_i}.
\]

Their Eq. 45 combines these segment integrals with each segment enthalpy rise
to calculate the path efficiency. This is precisely why a single ratio of
endpoint differences is not a general real-gas polytropic calculation: the
curved path integral must be resolved. The article recommends checking the
segment count/convergence for the chosen path shape. TurboCompGPU currently
does not implement the Taher–Evans cubic-segment method; its new diagnostic
instead solves a differential constant-ηp path from the EOS derivatives. The
Eckardt output still includes the total-to-total isentropic efficiency defined
above, with ηp emitted as a separate field. The published method and its assumptions
are documented in the cited article; the formula here is a short equation-level
summary, not a reproduction of its tables or standard text.

The implemented shooting path uses

\[
\frac{dT}{dp}=
\frac{1/(\eta_p\rho)-(\partial h/\partial p)_T}
     {(\partial h/\partial T)_p},\qquad
H_p=\int_{p_1}^{p_2}\frac{dp}{\rho}.
\]

It integrates (T) and (H_p) together with fixed-step RK4 in log pressure,
then bisects on constant ηp until the integrated path matches the actual
outlet total temperature. The reported efficiency is recomputed as
(H_p/(h_{03}-h_{01})). The calculation assumes a single-phase EOS path with
valid analytic derivatives and a constant polytropic efficiency over the path.
It is a diagnostic calculation and does not change the stage solution.

Complete bibliographic records and access/licensing notes are in
[`../references/README.md`](../references/README.md). The standards and books
are citation-only; no copyrighted text or PDFs were added.

## What is established versus still open

### Established by the current evidence

- The stage solver has equation-level conservation checks and explicit
  energy-residual reporting.
- The Eckardt benchmark's efficiency output is total-to-total and uses the
  selected thermodynamic model for the isentropic reference state.
- Published loss bundles are represented as diagnostic variants with their
  known geometry and closure proxies documented.
- Public papers provide useful evidence for mixture properties, compressor
  maps, loss-correlation behavior, and diffuser pressure measurements, but
  those sources do not all share the same machine, stations, gas, or output
  definitions.

### Required before claiming method/model validation

1. Reproduce a published real-gas polytropic numerical example, with
   source conventions, units, EOS, gas composition, and segment count recorded.
   The current constant-ηp integration already passes ideal-gas and
   pressure-dependent nonideal-gas analytic checks; do not call it PTC 10
   compliant until it is compared against the controlled code method.
2. Add an uncertainty ledger for each validation case: instrument/test
   uncertainty, property-model uncertainty, geometry/input uncertainty, and
   digitization uncertainty. Propagate them to pressure ratio, work, and each
   efficiency definition.
3. Score Avramenko diffuser data only after the model and data use matching
   diffuser geometry, station-average definitions, and static-pressure
   quantities. Keep unsteady/stalled regimes separate from steady-model
   qualification.
4. Keep the natural-gas compressor gate blocked until public matched geometry,
   gas composition, station definitions, operating data, and uncertainties are
   available. A property grid or synthetic Eckardt geometry is not a machine
   validation.

## Current acceptance statement

This repository has a stronger and more traceable reference basis and a checked
EOS-based constant-polytropic-efficiency diagnostic, but **method/model
validation remains incomplete**. The Eckardt
air-map residuals remain descriptive because digitization and test uncertainty
are not propagated; the PR mixture comparison fails its present GERG screen;
and no direct natural-gas stage map is available. A method citation, solver
convergence, or a close fit to the same data used to explore correlations is
not a validation pass.
