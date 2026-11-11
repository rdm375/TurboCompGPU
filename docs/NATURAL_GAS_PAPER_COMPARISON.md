# Natural-gas compressor paper comparison

This note compares the public evidence in Chen et al. (2025) and Lin et al.
(2026) with the current TurboCompGPU evidence recorded for remote `main` at
`2f8c993062b109ccbf652725befea476641c59d4`. It is an evidence comparison,
not a qualification claim.

## Published results and how they relate to this repository

| Source | Machine and method | Public result | Relation to TurboCompGPU |
|---|---|---|---|
| Chen et al. (2025), Tables 2-4 | A gas-storage centrifugal compressor; similarity-based performance-curve conversion; five field operating points and gas composition | Their converted values differ from reported field values by 0.32-0.57% in pressure ratio and 0.96-2.3% in power. Mean of the five rounded absolute errors is 0.386% for pressure ratio and 1.592% for power. | Useful as a natural-gas system/map benchmark and for a future mixture-EOS regression. It does not publish matched stage geometry or uncertainty and does not test a geometry-resolved meanline model. We have not run TurboCompGPU on this machine. |
| Lin et al. (2026), Sections 2-4 | A three-stage pipeline compressor; full-machine CFD baseline and factory-data residual correction | PR EOS vs GERG-2008 maximum deviations: Z 1.09%, density 1.11%, enthalpy 1.35%, heat-capacity ratio 1.72%. Corrected model maximum errors on the independent factory subset: outlet pressure 0.451%, temperature rise 0.918%, power 0.611%, polytropic head 0.481%, polytropic efficiency 0.899%. Post-overhaul field maxima: pressure 1.650%, power 3.048%, efficiency 1.708%. | Strong methodological reference for real-gas property checks, polytropic-head integration, independent evaluation, and separating fluid power from shaft power. Its pointwise measurements and reconstructed geometry are confidential, so its scores cannot be independently reproduced here. |
| TurboCompGPU at `2f8c993` | Eckardt A air map; PR pseudo-air and ideal-gas runs | PR-air pressure-ratio MAPE 1.328% (max absolute 3.221%); total-to-total efficiency MAE 4.766 percentage points (max absolute 7.663 pp). The ideal-gas results are nearly the same. | These are air-map results on a different single-stage machine. They are not directly comparable to either natural-gas study and do not demonstrate natural-gas accuracy. |

The Chen paper reports the actual and converted rows in Tables 2-3. Its
abstract describes conversion error as within 2%, but Table 3 reports a 2.3%
power error at operating condition 2; the comparison above preserves the
rounded table value and flags this source-level discrepancy. The raw
transcription is in
[`data/validation/chen_etal_2025_gas_storage_points.csv`](../data/validation/chen_etal_2025_gas_storage_points.csv).
The accompanying composition is in
[`data/validation/chen_etal_2025_gas_composition.csv`](../data/validation/chen_etal_2025_gas_composition.csv).
The source reports 0.24 mol% as `C5+` and does not resolve butane isomers;
these entries must remain explicit in any property-model setup. The five
published volume-flow values are retained in the paper's unit of 10^3 m^3/h;
the reference volume conditions are not specified in the tabulated data.

## What the current results establish

At the pinned source revision, the completed compressor comparisons use air
benchmark data. The PR-air property case uses a pseudo-pure dry-air model with
constant heat capacity. It is not the Chen or Lin natural-gas mixture. The
near-equality of the ideal-air and PR-air Eckardt predictions follows from
those near-ideal air conditions and the model setup; it says nothing about
compressibility corrections for methane-rich gas at 5.9-8 MPa.

ThermoGPU now provides a fixed-composition Peng-Robinson mixture model through
the same `thermomodel::Model` and `DifferentiableModel` interfaces used by the
stage solver. TurboCompGPU exposes the Lin et al. (2026) eight-species natural-gas
mixture as `pr-natural-gas-lin2026`. The component mole fractions and critical
properties are explicit; binary interaction coefficients are zero because the
paper does not publish its matrix, and ideal-gas component heat capacities are
constant screening values. The new property path therefore accepts natural-gas
composition in stage state calls and analytic property derivatives.

This closes the software capability gap, but does **not** create an aerodynamic
paper score. Lin et al. publish only aggregate PR-versus-GERG deviations and
withhold pointwise compressor data and geometry. Chen et al. provide field rows
and a composition, but not matched stage geometry or a defined reference-volume
basis for all flow values; additionally, its `C5+` lump and unresolved butane
isomers are not silently mapped to the Lin composition. A direct TurboCompGPU
compressor residual remains unavailable. The paper values are not calibration
targets, and no paper values were fitted.
The reproducible property-grid runner is `turbocompgpu_natural_gas_property_grid`.
It evaluates the 20 inlet combinations listed by Lin et al. (2026): five
pressures from 6.5 to 8.0 MPa and four temperatures from 20 to 35 °C. With the
published eight-species composition and the screening caloric assumptions
above, the PR model gives Z=0.8185–0.8721, density=51.01–70.31 kg/m³,
cp/cv=1.533–1.667, and cp=2.555–2.808 kJ/(kg K). The computed enthalpy range
(−107.5 to −50.5 kJ/kg) uses the model's explicitly declared zero ideal-gas
enthalpy at 298.15 K; it cannot be compared to absolute enthalpy without Lin's
reference-state convention.

Lin et al. report only aggregate PR-versus-GERG maximum/average deviations
(1.09/0.95% for Z, 1.11/0.96% for density, 1.35/0.88% for enthalpy, and
1.72/1.55% for specific heat ratio). They do not publish the 20 paired PR and
GERG values, so those paper scores cannot be independently reproduced. We
added a separate pointwise comparison against NIST's GERG-2008 implementation
using the same 20 T,p points and the Lin composition. With enthalpy shifted to a
common 298.15 K, 101.325 kPa state and specific heat ratio defined consistently
as cp/cv, the present PR setup gives maximum deviations of 2.914% for Z, 3.004%
for density, 19.512% for relative enthalpy, and 2.251% for cp/cv. The Z,
density, and cp/cv deviations exceed the corresponding Lin maxima. The
enthalpy percentage is not directly comparable because Lin et al. do not state
an enthalpy reference datum, while our comparison shifts each EOS to a shared
reference state. On a reference-invariant diagnostic, the PR-versus-GERG
discrepancy in the isothermal 6.5-to-8.0 MPa enthalpy change is 1.00–1.06
kJ/kg (5.9–7.1% of the GERG change) across 293.15–308.15 K. This is not a
compressor-work prediction. This is a failed EOS qualification screen, not a
reproduced paper score: the current PR setup assumes zero binary
interaction coefficients and constant ideal-gas component cp values, neither
of which is specified by the paper. The exact comparison and pointwise results
are recorded in
[`benchmarks/results-natural-gas/lin2026_pr_vs_gerg2008_summary.md`](../benchmarks/results-natural-gas/lin2026_pr_vs_gerg2008_summary.md).
The paper's “specific heat ratio” is compared as cp/cv; GERG's isentropic
exponent is a distinct property and is not substituted for cp/cv. The
one-component methane mixture limit matches the existing pure-fluid PR
backend to relative error below 2.3e-16 for Z, density, enthalpy, entropy,
and cp. At 303.15 K and 7.5 MPa, all six returned T,p derivatives agree
with centered finite differences to relative error below 2.1e-9; a
density-to-pressure round trip differs by 1.9e-9 Pa. The
`pr-natural-gas-lin2026` option on the Eckardt sweep selects the same
composition through the stage solver interface, but that sweep retains Eckardt
A geometry and its standard inlet boundary. It is an integration smoke path
only; its compressor residuals must not be interpreted as natural-gas
validation. A compiled run completed synthetic 14,000 rpm stage points from
0.25 to 2.00 kg/s. Reusing air-map operating flows caused an existing stage
physical check to reject a point, confirming the Eckardt air map is not a
natural-gas operating map.

## Reproducible comparison boundary

Chen et al. provide a small tabulated field dataset but not the detailed
machine geometry, station conventions, uncertainty, or reference-volume basis
needed to instantiate this solver and score a like-for-like prediction. Lin et
al. provide a detailed methodology and aggregate scores, but explicitly state
that the data are unavailable due to confidentiality. Both are still valuable:
Chen is a future system-level mixture/map check; Lin supplies property and
thermodynamic-definition checks and realistic expectations for a staged
validation workflow.

A valid direct comparison requires the same compressor, composition, operating
points, flow basis, station definitions, and reported output definitions. The
current PR mixture fails the external GERG-2008 property screen; resolving its
binary interaction and caloric inputs or replacing it with a validated
reference-quality mixture backend is required. The compressor gate also
remains blocked on an openly reusable natural-gas dataset with stage geometry
and matched measurements. Accelerator agreement is a separate numerical gate
and does not repair either physics gap.

## Sources

- Chen, H. et al. “The Influence of Variable Operating Conditions and
  Components on the Performance of Centrifugal Compressors in Natural Gas
  Storage Reservoirs.” *Energies* 18(15), 3930 (2025).
  [DOI](https://doi.org/10.3390/en18153930). Tables 2-4 provide the operating
  points, conversion results, and gas composition. Article data availability
  statement: the data supporting the findings are within the article.
- Lin, Y. et al. “Healthy-State Performance Modeling of a Multistage Natural
  Gas Centrifugal Compressor Using a CFD-Generated Baseline and Factory-Data
  Correction.” *Processes* 14(13), 2154 (2026).
  [DOI](https://doi.org/10.3390/pr14132154). The paper reports a PR/GERG-2008
  property comparison and aggregate factory/field errors; its data availability
  statement says the data are not public due to confidentiality restrictions.

Both publisher articles are open access under CC BY 4.0. Publisher PDFs are
archived in [`docs/references/`](references/README.md) with attribution and
checksums.
