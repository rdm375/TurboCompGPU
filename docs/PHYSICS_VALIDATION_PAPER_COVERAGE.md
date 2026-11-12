# Physics validation paper coverage

Updated: 2026-11-11

## Coverage estimate

The repository now has relevant literature for **6 of 7 physics-validation
domains (86%)**. This is a literature-coverage count: a domain counts when at
least one source provides a relevant equation, measured case, property
reference, or validation method. It does **not** mean 86% of the model is
validated or that the tests pass.

Of the same seven domains, **4 of 7 (57%) have accessible numerical reference
quantities that can be used for a reproducible comparison**. One of those,
the natural-gas mixture property comparison, is already run and currently
fails the reported Lin bounds. Another, the Avramenko diffuser dataset, is
transcribed but not yet compared. One domain—direct natural-gas full-stage
prediction against a geometry-matched public test—has no usable paper dataset.

By the stricter measure of comparisons already run against a relevant
reference, the repository has **2 of 7 domains (29%) evaluated**: a natural-gas
property screen and the Eckardt A air-map screen. The property screen fails
published bounds; the air-map screen still has unresolved efficiency residuals.
Thus, “evaluated” does not mean “passed.”

Accelerator CPU/GPU/TPU equivalence is not counted as a paper-coverage domain;
it is a software parity experiment against a frozen CPU reference.

## Domain inventory

| Physics domain | Paper coverage | Numerical reference availability | Present status |
|---|---|---|---|
| Mixture EOS and real-gas properties | Lin et al. (2026), Srinivasan et al. (2024), plus GERG-2008 reference | Yes; Lin's exact paired EOS values are not public, but our own GERG comparison uses the same published mixture and state grid | Implemented PR screen fails reported Z, density, and cp/cv maxima; enthalpy datum is not comparable |
| Thermodynamic work and efficiency definitions | Srinivasan et al. (2024), Lin et al. (2026), Chen et al. (2025) | Yes for published PTC outputs and selected field points; conventions must be matched | Definitions are documented; a matching PTC-10 implementation/reference comparison is still needed |
| Impeller velocity triangles and empirical loss correlations | Kovář et al. (2021), Yang et al. (2023), Eckardt and related published loss sources | Partial: correlations and map plots/tables exist, but some results are digitized or coefficient-optimized | Useful implementation cross-checks; fitted coefficients are not independent prediction evidence |
| Diffuser pressure recovery and off-design behavior | Avramenko et al. (2015), Jaatinen-Värri et al. (2014), Dewar et al. (2019), Rossbach et al. (2015) | Yes: Avramenko publishes eight station regimes; Dewar has nine map points and tap locations; Rossbach Fig. 9 has measured component recovery curves | Avramenko and Rossbach inputs are transcribed at different precision levels; no direct score with the current Eckardt-only executable |
| Air compressor stage/map performance | Eckardt A, Choi–Kang, Kovář et al. (2021/2023), Dewar et al. (2019), Rossbach et al. (2015), Methel (2016) | Yes, with tabulated and figure-digitized sources on different machines | Screened only on Eckardt A; its efficiency-slope and magnitude residuals remain material. Other source geometries need case adapters |
| Measurement uncertainty and test interpretation | Srinivasan et al. (2024), Chen et al. (2025), Avramenko et al. (2015) | Partial: selected measurement practices and points are published, but uncertainty is not consistently sufficient for our acceptance limits | Our own acceptance criteria and uncertainty propagation still need completion |
| Natural-gas compressor performance with matched stage geometry | Chen et al. (2025), Lin et al. (2026) discuss real machines, but do not provide the matched public geometry and pointwise stage data required here | No usable blind stage-map dataset | **Uncovered; blocks natural-gas aerodynamic validation** |

The 6/7 count includes literature context for the natural-gas stage domain only
if it supplies a usable direct benchmark. Chen and Lin do not, so this domain
is the single uncovered item. Their papers remain useful as measurement and
thermodynamic-method references.

## What each newly added paper contributes

- **Srinivasan, Kurz, and David (2024), GT2024-125064:** natural-gas
  composition, field checkpoints, EOS comparison, and PTC-10 polytropic
  calculations. Use it to validate mixture thermodynamic and efficiency
  calculations, keeping polytropic efficiency separate from stage total-to-total
  efficiency. It is a multistage machine and is not a matched stage loss-map
  case. The paper excludes GERG and AGA from its EOS comparison.
- **Avramenko, Frolov, and Hämäläinen (2015):** eight experimental regimes for a
  low-pressure single-stage compressor with a vaneless diffuser, including
  station pressure measurements. Ratios are section-average absolute static
  pressures. Use only for the appropriate component-pressure definition and
  high-flow regime; low-flow regimes show unsteady stall behavior.
- **Yang, Liu, and Zhao (2023):** 1D off-design loss-model treatment with
  splitter-blade cases, including Eckardt-O and SRV2-O. Use it to compare
  correlation implementation and trend behavior. Their optimized-coefficient
  results are calibration, not an independent accuracy target.

## References and artifacts

Paper citations, access links, license notes, and the archived Yang et al. PDF
are in [`docs/references/README.md`](references/README.md). Avramenko's eight
regimes are in
[`data/validation/avramenko_etal_2015_diffuser.csv`](../data/validation/avramenko_etal_2015_diffuser.csv).

The remaining physics-validation gates are recorded in
[`PHYSICS_VALIDATION_EXIT_GATE.md`](PHYSICS_VALIDATION_EXIT_GATE.md). Paper
coverage does not change those gate statuses.

## Additional non-NASA datasets found (2026-11-11)

Dewar et al. (2019) supplies nine operating points and 12 diffuser tap
positions; the tables are transcribed in `data/validation/dewar_etal_2019_*`,
but Figure 8 still needs calibrated digitization. Rossbach et al. (2015)
supplies full-stage and component-level air-rig measurements in figures. Its
Figure 9 experimental pressure-recovery curves are approximately digitized in
`data/validation/rossbach_etal_2015_diffusion_pressure_recovery_digitized.csv`;
the trend screen passes, but this is not a TurboCompGPU run. Methel (2016)
supports qualitative diffuser-design and clearance conclusions; the thesis
download was blocked in this environment, so it is not yet transcribed.
DATED supplies outputs from a 1D meanline model and can support cross-model
checks only after a matching geometry and license are established. See
`docs/EXTERNAL_DATASETS_NON_NASA.md` for source and reuse details.
