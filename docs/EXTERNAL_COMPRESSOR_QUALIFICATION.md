# External compressor qualification: Eckardt B, Schiffmann–Favrat, RADIVER, Sandia

## Eckardt B: structured speed-line holdout

The Eckardt B benchmark reuses the existing map evaluator with the B rotor
geometry from Meroni et al., Table 5, and the 17 tabulated pressure-ratio points
from Bayomi et al., Table 5. The geometry is distinct from Eckardt A, but both
rotors belong to the same Eckardt experimental family. The inlet annulus is
represented by a circular-equivalent meanline area, and the solver has no
Eckardt B-specific inlet incidence map. This is an external map screen, not an
independent rig or a full geometry-resolved qualification.

The plain map run compares measured pressure ratio against the geometry-based
baseline. It remains a diagnostic screen: its default slip correlation is not
an accuracy pass/fail test.

The default correlation run converges at all 17 points but has 10.345% MAPE
and 20.576% maximum error. A separate qualification runner fits an
**empirical effective slip/head parameter** on the 10,000, 12,000, and
14,000 rpm lines, then scores the untouched 16,000 rpm line. It does not call
the fitted parameter a physical slip measurement. The runner explicitly
selects Eckardt B geometry for every fit and prediction. Recomputed against the
pushed B geometry, the fit is 0.74572749; training MAPE is 0.921%, and the
held-out 16,000 rpm line has 3.785% MAPE and 5.151% maximum absolute error.

The predeclared default acceptance gate is all holdout points converged and
maximum absolute pressure-ratio error no greater than 7%. The 7% value is a
provisional project tolerance informed by the maximum B pressure-ratio
deviation reported for a separate published mean-line model. It is not a
measurement uncertainty bound. This PASS qualifies prediction of one B speed
line after same-rotor calibration; it does not qualify an unfitted model or an
independent compressor family.

Run the baseline screen and held-out qualification through CTest:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_eckardt_impeller_a \
           turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^(eckardt_[ab]_external_map_screen|eckardt_b_heldout_speedline_qualification)$'
```

The baseline CSV is `build-eckardt-validation/eckardt_b_qualification.csv`.
The fitted held-out CSV and report are
`build-eckardt-validation/eckardt_b_speedline_holdout.csv` and `.md`.

The station-matched Choi–Kang Figure 2 comparison is also runnable as a
three-point Eckardt A screen:

```sh
cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8
ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^eckardt_a_choi_kang_station_matched_screen$'
```

That CTest verifies the comparison runner completes; its approximate
digitization and efficiency results are reported as a screen, not an accuracy
gate.

## Schiffmann–Favrat: R134a design-point smoke screen

Figure 10 in Schiffmann and Favrat (2009) contains measured speed lines from
150 to 210 krpm and internal isentropic efficiency contours for an R134a
compressor. The figure-digitization template remains empty at
`data/validation/schiffmann_favrat_fig10_template.csv`; it has no invented
coordinates.

There is now a runnable one-point screen using the summary point in Meroni's
thesis (180 krpm, 0.039 kg/s, PR 2.35, efficiency 0.80) and the geometry in
its Appendix C. This is **not** a digitization of Figure 10, map qualification,
or acceptance test. R134a is parameterized in generic Peng–Robinson with a
constant-Cp caloric approximation; that property setup is not reference-
qualified. The stage smoke model supplies zero internal loss, so its
efficiency is an idealized 1.0 and cannot be compared as a compressor
efficiency prediction. Its current PR mismatch is also large. Use this runner
only to verify EOS, geometry, and solver wiring; do not tune against this one
secondary-source summary point.

There is a source-condition discrepancy to resolve before a performance
comparison: the summary-point CSV uses 265 K, while a later published
validation table lists 283.15 K for the EPFL compressor. Confirm the original
run conditions and station definitions, then digitize Figure 10 with
uncertainty and use a reference-quality R134a property model and nonzero loss
inputs.

Run it with:

```sh
cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_schiffmann_favrat_design_point --parallel 8
ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^schiffmann_favrat_design_point_smoke$'
```

## RADIVER: local licensed-map runner for the vaneless case

RADIVER has multiple tested diffuser configurations, including vaneless and
vaned variants. Choose the vaneless-diffuser configuration to match the
current model; the vaned variant is outside the current model. The repository
includes an empty CSV schema at
`data/validation/radiver_map_template.csv` and a local-map runner. Do not commit
licensed measurements or digitized derivatives. Keep the authorized map
outside the repository and record its source run, measurement plane, operating
condition, and uncertainty in a local manifest. The runner accepts explicit
`impeller_exit_total` and `vaneless_exit_total` station labels; it rejects
other labels rather than comparing unmatched quantities.

RWTH describes RadiVer 2 as available for CFD validation, with a signed
agreement that forbids sharing access to or distributing the supplied test
case. Request the data from RWTH and keep the received files and derivatives
out of the public repository unless the agreement explicitly allows them.
Configure a local data path, build, and run the screen as follows:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu" \
  -DTURBOCOMPGPU_RADIVER_LOCAL_MAP=/path/to/authorized/radiver_map.csv

cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_radiver_local_map --parallel 8
ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^radiver_authorized_local_map_screen$'
```

The runner uses public meanline geometry summaries plus provisional blockage
and slip assumptions. It is a solver screen with no accuracy gate; keep the
simplifications visible alongside the residuals.

## Sandia sCO2: runnable impeller-only screen

ThermoGPU's Peng–Robinson class is a generic pure-fluid implementation, so it
can be parameterized for CO2. The Sandia check constructs pure CO2 using
critical constants and a NIST Shomate ideal-gas heat-capacity correlation. This
is API support, not a validated CO2 property package: Peng–Robinson differs
substantially from Span–Wagner near the pseudo-critical region.

The source table's inlet pressure values (76.76, 79.20, etc.) are in **bar**,
despite a secondary paper labeling the column MPa. Monje's Table 10 explicitly
identifies bar, and Sandia's report gives the inlet near 7.687 MPa. The CSV
stores the converted SI pressures (7.676 MPa, etc.) to prevent a factor-of-ten
mistake. It also records the impeller-exit static pressure P3 and the vaned
diffuser outlet P05/T05 values reported by Monje. Only P3 is compared here;
P05/T05 belong to a vaned diffuser that TurboCompGPU does not model.

Run the component screen:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_sandia_sco2_impeller --parallel 8

ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^sandia_sco2_impeller_component_screen$'
```

The runner converts inlet total state to static state using the PR model,
published inlet annulus geometry, and continuity. It then compares the model's
impeller-exit static pressure with the source's P3. The benchmark uses the
published impeller dimensions, an estimated blade blockage, and a provisional
slip factor; it does not model the vaned diffuser. Therefore this is a
component screen, not a qualification or a like-for-like reconstruction of
Monje's meanline model. There is no accuracy gate. On the current run, all 12
operating points converged; P3 MAPE was 4.92%, with 9.68% maximum absolute
error. These errors include property, inlet-geometry, blockage, slip, and
impeller-model differences and should not be assigned to one cause.

The executable retains two spot checks against Span–Wagner reference tables at
8 MPa:

| State | Density error | Cp error |
|---|---:|---:|
| 305 K | -12.33% | +20.56% |
| 310 K | +1.07% | -12.70% |

For a broader diagnostic, install CoolProp into the Python environment selected
by CMake. Its CO2 HEOS backend uses the Span–Wagner formulation. The property
comparison exports Peng–Robinson states at all 12 Sandia inlet conditions and
over a nearby temperature-pressure grid. It reports separate errors for the 12
Sandia operating states, the diagnostic grid, and the enthalpy/entropy datum
anchor. Only the operating-state group describes the tabulated compressor map;
the nearby grid is a diagnostic stress sample. There is no property accuracy
gate.

```sh
python3 -m pip install CoolProp
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"
cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_sandia_sco2_impeller --parallel 8
ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^sandia_co2_properties_against_heos$'
```

The same CoolProp installation can also run the compressor screen with HEOS
inside the C++ solver. This is an optional executable: CMake must find the
selected Python interpreter, its embedding development files, and CoolProp in
that interpreter. Reconfigure after installing CoolProp so CMake can detect
those dependencies. If available, run:

```sh
cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_sandia_sco2_impeller_heos --parallel 8
ctest --test-dir build-eckardt-validation --output-on-failure \
  -R '^sandia_sco2_impeller_heos_component_screen$'
```

The HEOS and Peng–Robinson component screens use the same measured points,
geometry, inlet reconstruction, and provisional slip factor. The HEOS run has
no acceptance threshold; it checks convergence and reports the outlet-pressure
residuals when the EOS is changed. It does not qualify the geometry or loss
model, and its property reference convention differs from the PR model.

An exploratory local run produced valid impeller results at 9 of 12 points.
On those same nine points, HEOS P3 MAPE was 8.10% (maximum 16.09%), versus
4.56% (maximum 9.68%) for PR. The three HEOS inlet reconstructions that did
not complete were 29,888, 55,000, and 60,000 rpm; the entropy-constrained
static-inlet solve could not bracket enthalpy with the current geometry and
inlet reconstruction. This indicates that replacing PR alone does not close
the Sandia residuals; the comparison still combines property, inlet-state,
slip, and impeller-model effects. Treat these numbers as a diagnostic run, not
a qualified EOS or compressor result.

In the grouped property audit, PR density and Cp MAPE at the 12 Sandia inlet
states were 12.95% and 11.76%, respectively (maximum absolute errors 14.14%
and 37.84%). The separate 30-state diagnostic grid gave density and Cp MAPE of
8.68% and 19.84%. These groups should remain separate: only the first uses the
reported map inlet conditions.


To test whether one empirical effective slip/head scalar explains the P3
residual trend, the HEOS executable has an optional bounded calibration mode.
It fits on the five points at or below 39,000 rpm and evaluates the seven
higher-RPM points as holdout data. Since the table has at most one point at
each RPM, this is a low-speed/high-speed split, not a holdout of complete
measured speed lines. The inlet static states are computed once and reused
during fitting. The objective is training-set mean squared relative P3 error;
the default scalar bounds are 0.60–0.99. The output compares the baseline 0.84
and fitted scalar on both sets without a pass/fail threshold. This parameter is
an empirical head correction, not a measured blade slip factor.

The bounded nonlinear fit requires repeated operating-point solves and can
take several minutes on the CPU. Run it after building the optional HEOS
executable:

\`\`\`sh
./build-eckardt-validation/turbocompgpu_validate_sandia_sco2_impeller_heos \\
  --fit-effective-slip \\
  data/validation/sandia_sco2_report_table.csv \\
  39000 \\
  build-eckardt-validation/sandia_sco2_heos_effective_slip_fit.csv
\`\`\`

The CSV contains baseline and fitted errors for every row, tagged as training
or holdout. This experiment tests transfer of the low-speed calibration to
higher speeds; it still does not qualify the estimated geometry or impeller
loss model.

These comparisons do not make the impeller screen a CO2 qualification. Do not
tune the air loss model against the Sandia stage outlet table. A stronger
near-critical CO2 claim requires a validated Span–Wagner-class model in the
compressor calculations, not only a reference-comparison script.

## Evidence summary

| Dataset | Source data in repo | Model comparison now | Main limit |
|---|---|---|---|
| Eckardt B | 17 tabulated points | Yes: baseline screen plus 16 krpm held-out line gate | Same Eckardt family; fitted effective parameter is rotor-specific calibration |
| Choi–Kang / Eckardt A | 3 approximate digitized points | Yes: station-matched CTest screen | Approximate digitization; one speed only |
| Schiffmann–Favrat | One secondary-source design summary plus empty Figure 10 template | Fresh smoke run converged but predicted PR 3.0732 vs 2.35 (+30.78%) and eta=1.000 vs 0.800 (+20.0 pp) | Poor agreement; constant-cp PR properties, reconstructed geometry, and provisional slip make this a smoke test, not qualification |
| RADIVER | Empty licensed-data schema | Local-map runner available; requires authorized input | Public meanline geometry and provisional slip/blockage; no committed test points |
| Sandia sCO2 | 12 map points plus P3/P05/T05 references | Fresh PR/Shomate P3 run: 12/12 converge, 4.917% MAPE, 9.684% max; optional HEOS data remain historical and unavailable in this runtime | PR spot Cp errors reach 20.56%; estimated geometry/blockage/slip; vaned outlet not modeled; no full-stage CO2 qualification |

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
- The secondary summary point is transcribed from Meroni et al. (2018),
  *Design of centrifugal compressors for heat pump systems*, Table 4.4 and
  Appendix C. The inlet-temperature discrepancy above remains unresolved.
- A later EPFL comparison reports the compressor inlet at 283.15 K: Meroni
  et al. (2022), “The Effect of Size and Working Fluid on the Multi-Objective
  Design of High-Speed Centrifugal Compressors,” DOI: 10.1016/j.ijrefrig.2022.06.023.
- RADIVER license agreement from RWTH IST:
  <https://www.ist.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaoncimv>.
  It permits the licensee's own R&D use and prohibits distributing case
  contents to third parties; follow the signed agreement for the copy received.
- Jung et al. (2012), “Effects of a Nonuniform Tip Clearance Profile on the
  Performance and Flow Field in a Centrifugal Compressor,” Table 1 and 80%
  speedline discussion, DOI: 10.1155/2012/340439. Geometry and inlet conditions
  for the local RADIVER runner.
- Wright, S. A. et al. (2010), *Operation and Analysis of a Supercritical CO2
  Brayton Cycle*, SAND2010-0171, DOI: 10.2172/984129; Sandia measurement maps.
- Sandia point values reproduced in “A New Method for Impeller Inlet Design of
  Supercritical CO2 Centrifugal Compressors in Brayton Cycles,” *Energies*
  13(19), 5049 (2020), Table 3, DOI: 10.3390/en13195049. The original Sandia
  report remains the primary source for the experiment.
- Monje, B. (2014), *Design of supercritical carbon dioxide centrifugal
  compressors*, Tables 9–10. Table 10 identifies the Sandia inlet-pressure
  column in bar, reports P3 and P05 in bar, and provides T05 for three cases.
- Span, R. and Wagner, W. (1996), “A New Equation of State for Carbon Dioxide
  Covering the Fluid Region from the Triple-Point Temperature to 1100 K at
  Pressures up to 800 MPa,” *J. Phys. Chem. Ref. Data*, 25, 1509–1596. The two
  8 MPa property spot checks use the 305 K and 310 K isobar entries.
- CoolProp 8 documentation for CO2 identifies the HEOS backend's Span–Wagner
  equation of state: <https://coolprop.org/fluid_properties/fluids/CarbonDioxide.html>.
