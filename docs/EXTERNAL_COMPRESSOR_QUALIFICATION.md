# External compressor qualification: Eckardt B, Schiffmann–Favrat, RADIVER, Sandia

## Implemented now: Eckardt B map screen

The Eckardt B benchmark reuses the existing map evaluator with the B rotor
geometry from Meroni et al., Table 5, and the 17 tabulated pressure-ratio points
from Bayomi et al., Table 5. The geometry is distinct from Eckardt A, but both
rotors belong to the same Eckardt experimental family. The inlet annulus is
represented by a circular-equivalent meanline area, and the solver has no
Eckardt B-specific inlet incidence map. This is an external map screen, not an
independent rig or a full geometry-resolved qualification.

The benchmark compares the measured full-stage pressure ratio against the
current lossless vaneless-diffuser model. It reports every point, solver
failures, bias, MAPE, and maximum absolute error. No coefficients are fitted to
these measurements.

Run the comparison before interpreting the model error. The test executable's
zero/nonzero exit status reports whether any points solved; it does not encode
a physical acceptance threshold. The map uses the same broader experimental
family as Eckardt A, so it is an additional rotor-map screen rather than a
fully independent rig qualification.

Run it through CTest:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_eckardt_impeller_a --parallel 8

ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^eckardt_b_external_map_screen$'
```

CTest writes pointwise predictions to
`build-eckardt-validation/eckardt_b_qualification.csv`.

## Schiffmann–Favrat: digitized map, evaluator blocked on fluid support

Figure 10 in Schiffmann and Favrat (2009) contains measured speed lines from
150 to 210 krpm and internal isentropic efficiency contours for an R134a
compressor. A transcription template is provided at
`data/validation/schiffmann_favrat_fig10_template.csv`; it intentionally has no
invented coordinates. Populate it only from a documented digitization, keeping
pixel coordinates and uncertainty in a companion manifest.

The source inlet is R134a at approximately 265 K and 1.65 bar. A meaningful
prediction requires a validated R134a property model, its source-matched
internal-efficiency definition, and the small compressor's measured impeller
and diffuser geometry. The current production validation executable is
configured only for air and cannot produce a qualified comparison here.

## RADIVER: licensed local dataset; select the diffuser configuration

RADIVER has multiple tested diffuser configurations, including vaneless and
vaned variants. Choose the vaneless-diffuser configuration to match the
current model; the vaned variant is outside the current model. The repository
includes only an empty CSV schema at
`data/validation/radiver_map_template.csv`. Do not commit licensed measurements
or digitized derivatives unless the applicable terms permit it. Keep any
authorized source files outside version control, and record source run,
measurement plane, operating condition, and uncertainty in a local manifest.

The RWTH license permits use for the licensee's own R&D and prohibits
distribution of the test-case contents to third parties. The first runnable
RADIVER test should consume a local, authorized export and leave raw and derived
licensed measurements outside the public repository. Its useful quantities are
the impeller-exit and vaneless-diffuser stations where the reported planes
match the implemented model.

## Sandia sCO2: source points recorded, full-stage comparison blocked

`data/validation/sandia_sco2_report_table.csv` records 12 experimental
operating points reproduced in a published comparison table. It includes
speed, inlet total state, mass flow, and measured total pressure ratio. The
Sandia stage uses a vaned diffuser; therefore those outlet ratios cannot be
compared as if they were vaneless-stage predictions. The near-critical CO2
inlet conditions also demand a CO2-calibrated real-gas model and correct
total-to-static conversion. The current generic Peng–Robinson configuration
and vaneless-stage model do not yet qualify this dataset.

Use the Sandia points to build two later checks: (1) a total-state and CO2 EOS
property regression around the measured inlet conditions, then (2) a
component-matched impeller-only comparison if an experimental impeller-exit
measurement can be sourced. Do not tune the current air loss model against
the vaned-stage outlet table.

## Evidence summary

| Dataset | Source data in repo | Model comparison now | Main limit |
|---|---|---|---|
| Eckardt B | 17 tabulated points | Yes, CTest map screen | Same family as Eckardt A; approximate inlet annulus |
| Schiffmann–Favrat | Figure digitization template | No | R134a properties, geometry, and station-matched evaluator |
| RADIVER | Empty private-data schema | No | Licensed files absent; select and match a diffuser configuration |
| Sandia sCO2 | 12 reported experimental points | No full-stage comparison | Near-critical CO2 EOS and vaned diffuser unsupported |

## Sources

- Bayomi, N. N., Abdel Hafiz, A., and Hegazi, A. A. (2013), “Centrifugal
  Compressor Map Prediction and Modification,” Table 5. The table reports
  measured Eckardt rotor maps and is used here for Eckardt B point values.
- Meroni, A. et al. (2018), “Design of centrifugal compressors for heat pump
  systems,” Table 5, DOI: 10.1016/j.apenergy.2018.09.210. Geometry and inlet
  conditions for Eckardt B and Schiffmann–Favrat.
- Schiffmann, J. and Favrat, D. (2009), “Experimental investigation of a
  direct driven radial compressor for domestic heat pumps,” Fig. 10,
  DOI: 10.1016/j.ijrefrig.2009.07.006.
- RADIVER source case description: Dambach et al. and related RWTH publications;
  obtain and follow the dataset's current license before loading or copying
  data.
- Wright, S. A. et al. (2010), *Operation and Analysis of a Supercritical CO2
  Brayton Cycle*, SAND2010-0171, DOI: 10.2172/984129; Sandia measurement maps.
- Sandia point values reproduced in “A New Method for Impeller Inlet Design of
  Supercritical CO2 Centrifugal Compressors in Brayton Cycles,” *Energies*
  13(19), 5049 (2020), Table 3, DOI: 10.3390/en13195049. The original Sandia
  report remains the primary source for the experiment.
