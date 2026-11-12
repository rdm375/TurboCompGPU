# Non-NASA external compressor datasets: acquisition and validation status

Reviewed 2026-11-11. NASA HECC is intentionally excluded from this campaign.

## What is included now

The repository includes the open-access Dewar et al. (2019) paper and a
source-labeled transcription of its two numerical tables that are directly
useful as validation inputs:

- `data/validation/dewar_etal_2019_diffuser_operating_points.csv`: nine
  measured operating conditions (speed, mass flow, total-to-total pressure
  ratio), plus source conditions and the reported design-point reference.
- `data/validation/dewar_etal_2019_diffuser_tap_locations.csv`: all 12
  diffuser radial pressure-tap positions and wall-side labels.
- `data/validation/rossbach_etal_2015_diffusion_pressure_recovery_digitized.csv`:
  approximate experimental points read from Figure 9, with a digitization
  uncertainty estimate and source figure recorded per row.
- `docs/references/rossbach-etal-2015-process-compressor-test-rig.pdf`: the
  proceedings paper, reproduced unchanged under Euroturbo's CC BY-NC-ND 4.0
  policy.

The paper reports that experimental static pressure was measured at those
12 positions for all nine conditions. Figure 8 plots the experimental and CFD
pressure ratios and reports their differences; these plot values have **not**
been transcribed into per-tap rows yet. Its plot-only data need to be
machine-digitized with a digitization uncertainty before use as exact scores.

The design-point result is `eta_tt=79.8%`, `PR_tt=2.36`, and `mdot=1.79 kg/s`;
the paper reports relative errors of 0.81%, 0.35%, and 0.54%, respectively.
The paper's nine-point map has no efficiency values for all nine points.

## Candidate-by-candidate disposition

| Source | Evidence and access | What we can validate | Current disposition |
|---|---|---|---|
| Dewar et al. (2019), LUT high-speed air compressor | Open-access paper archived below under CC BY. Tables 1-4 give geometry summaries, a design point, map operating points, and tap positions. Figure 8 gives plotted static-pressure profiles. | Overall pressure ratio at nine conditions; design-point efficiency; diffuser static-pressure radial trend. | Operating map and tap locations transcribed. Figure 8 pointwise measurements still need uncertainty-qualified digitization. No TurboCompGPU score: current external-map executable hard-codes Eckardt A geometry and inlet conditions. |
| Methel (2016), Purdue CSTAR diffuser comparison thesis | Open-access thesis landing page. It compares airfoil versus wedge vaned diffusers and open versus tight impeller clearances. | Diffuser pressure recovery/effectiveness, clearance trends, and overall map changes. | Purdue's download returned HTTP 403 in this retrieval environment; no thesis PDF or numerical values were obtained. The repository labels it open access but exposes no reuse license on the record; link to the source and check reuse terms before mirroring. The abstract supports qualitative conclusions only. |
| Rossbach et al. (2015), RWTH high-flow test rig | Euroturbo proceedings paper; working fluid is air despite “process gas” in the title. The paper is distributed under CC BY-NC-ND 4.0; unchanged PDF archived below. It contains a design-speed stage map and diffuser/return-channel pressure measurements. | Diffuser pressure recovery trend versus dimensionless flow (Fig. 9), radial diffuser pressure distribution at design point (Fig. 6), and full-stage map trends (Figs. 3-4). | Figure 9 experimental markers have been approximately digitized and checked against the paper's stated qualitative trends. The plot read-off is not raw test data. No TurboCompGPU score: the fixed Eckardt A stage cannot represent this rig; its bulk diffuser model also cannot be compared directly with hub/shroud wall taps without a measurement operator. |
| Picard, Schiffmann, and Ahmed (2023), DATED | Zenodo dataset, 22 million predicted cases from a 1D meanline model, including maps for 14 known geometries and train/test partitions. | Cross-model numerical regression and broad design-space/trend checks, after matched geometries and conventions are selected. | Not experimental evidence. The visible Zenodo license field is blank, so no data files are redistributed. A 9.3 MB `dated_real_map` subset was inspected locally for schema only; do not call DATED agreement physical validation. |

## What the current implementation prevents

The checked-in `turbocompgpu_eckardt_diffuser_friction_sweep` program uses a
fixed Eckardt A geometry and assumes an Eckardt-specific inlet annulus,
blade count, slip relation, and loss-geometry inputs. Running it against
Dewar, Purdue, or RWTH operating points would produce numbers, but they would
not be like-for-like predictions. This campaign therefore does not create or
publish such residuals.

To make these sources executable validation cases, the next engineering step
is a case-driven benchmark interface that receives (and records) geometry,
gas, station definitions, boundary states, measurement operators, and
loss-model choices per case. Then:

1. Digitize Dewar Figure 8's experimental markers at the 12 taps, retaining
   speed, flow, tap radius, wall side, source figure, and estimated digitization
   uncertainty for every row.
2. Add a Dewar case definition from the paper's geometry and test conditions;
   score total-to-total pressure ratio at all nine points, design-point
   efficiency, and tap static pressure separately.
3. Retrieve Methel's thesis from Purdue outside this environment, inspect its
   reuse terms, and extract numerical map and component pressure-recovery
   results while separating clearance variants and diffuser types.
4. Re-digitize Rossbach Figures 3, 4, and 6 with calibrated axes and
   uncertainty bounds; obtain a machine-readable test-case definition if the
   authors' announced open test case is available.
5. Use DATED only for implementation cross-checks after a geometry-matched
   subset, model conventions, and license are established.

No loss parameters should be tuned on the same operating points used for the
reported validation score. Use one speed line or condition subset for
calibration only if the result is explicitly labeled calibrated, and retain a
separate holdout for predictive validation.

## References

- Dewar, B., Tiainen, J., Jaatinen-Värri, A., Creamer, M., Dotcheva, M.,
  Radulovic, J., and Buick, J. M. (2019). “CFD Modelling of a Centrifugal
  Compressor with Experimental Validation through Radial Diffuser Static
  Pressure Measurement.” *International Journal of Rotating Machinery*, 2019,
  Article 7415263. DOI: [10.1155/2019/7415263](https://doi.org/10.1155/2019/7415263).
  CC BY, publisher PDF archived with attribution.
- Methel, C.-T. J. (2016). “An experimental comparison of diffuser designs in a
  centrifugal compressor.” Purdue University MS thesis.
  [Repository record](https://docs.lib.purdue.edu/open_access_theses/796/).
- Rossbach, T., Rube, C., Wedeking, M., Franz, H., and Jeschke, P. (2015).
  “Performance Measurements of a Full-Stage Centrifugal Process Gas Compressor
  Test Rig.” Euroturbo ETC2015-084. The test medium is air.
  [Paper record and download](https://www.euroturbo.eu/publications/proceedings-papers/etc2015-084/).
  The paper is licensed CC BY-NC-ND 4.0. Archived PDF SHA-256:
  `16e54539cc1d1b1599efc5fe786d9337efa0cca37b57f255b9d537acad313a15`.
- Picard, C., Schiffmann, J., and Ahmed, F. (2023). “DATED: A Dataset of
  Centrifugal Compressors.” Zenodo v1.0.0.
  [DOI 10.5281/zenodo.8200792](https://doi.org/10.5281/zenodo.8200792).
