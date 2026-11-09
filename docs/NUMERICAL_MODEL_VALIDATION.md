# Numerical compressor model validation campaign

Reviewed: 2026-11-09
Repository baseline: 3980363d090d3b532cf50d437f8b90cfc9a0fa21

## Scope and claim levels

This campaign validates the numerical compressor model component by component;
it does not qualify the ReLU surrogate. Keep four claims distinct:

1. **Numerical verification:** equations, derivatives, conservation,
   independent formulations, tolerances, and solver behavior.
2. **Component validation:** agreement with measurements of the same quantity
   at matched conditions and measurement planes.
3. **Calibrated prediction:** fitting on declared training data and scoring on
   untouched holdout data.
4. **System integration:** component/network closure and transient behavior.
   Passing an integration test does not validate compressor physics against a
   rig.

## Coverage map

| Area | Existing numerical evidence | Available measurement comparison | Limit |
|---|---|---|---|
| EOS and thermodynamic states | compressor_ideal_gas, compressor_peng_robinson, peng_robinson_derivatives, station_thermodynamics, isentropic_sound_speed | Sandia CO2 state comparison to CoolProp HEOS/Span–Wagner when dependencies are installed | PR has material near-critical CO2 property error; property agreement is not compressor validation |
| Nonlinear solver and operating point | nonlinear_residual, nonlinear_solver, nonlinear_solver_extended, nonlinear_solver_failure_paths, operating_point, compressor_failures | Eckardt map runners | Convergence does not establish physical admissibility |
| Impeller | impeller_aerodynamics, impeller_prewhirl, impeller_continuity, impeller_continuity_independent, impeller_conservation, impeller_system, impeller_system_loss, incidence_loss | Eckardt A/B maps; Sandia impeller-exit static pressure P3 | Inlet geometry, slip, blockage, loss closures, and station definitions affect residuals |
| Diffuser and choking | vaneless_diffuser, diffuser_characteristic, vaneless_diffuser_friction_verification, vaneless_diffuser_friction_analytical, vaneless_diffuser_variable_width_verification, diffuser_wall_geometry, choking_independent_validation, geometry_coupled_choking | Choi–Kang station-matched points; optional authorized RADIVER local map | No broad same-machine diffuser qualification |
| Complete stage | complete_stage_reference, complete_stage_friction_verification, complete_stage_recovery_geometry, complete_stage_reduced_recovery, impeller_diffuser_integration | Eckardt A/B and Choi–Kang screens | No geometry-complete independent full-stage dataset currently establishes broad accuracy |
| Maps and control interfaces | calibrated_map_interpolation, pipeline_compressor_control | TUM IGV summaries are cross-machine software/data-pipeline evidence | Selected-machine vane, valve/cooler, and surge maps are missing |
| Pipe coupling | Optional two_pipe_compressor_coupling | No matched measured pipe/compressor transient in the repository | Integration closure is numerical evidence only; actuator dynamics and thermal transport are open |

Test identifiers are from CMakeLists.txt at the baseline commit. Optional
dependencies and local licensed inputs affect which tests CMake registers.
passage_turning_sensitivity is intentionally disabled and must not be counted
as passed.

## Campaign order

### 1. Freeze the regression baseline

Use the same ThermoGPU install as prior Dell runs:

~~~sh
cmake -S . -B build-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-validation --parallel 8
ctest --test-dir build-validation --output-on-failure
~~~

Record commit SHA, compiler/version, configure options, optional dependencies,
enabled/disabled tests, and full output.

### 2. Run core numerical verification

~~~sh
ctest --test-dir build-validation --output-on-failure \
  -R '^(compressor_ideal_gas|compressor_peng_robinson|peng_robinson_derivatives|station_thermodynamics|isentropic_sound_speed|nonlinear_residual|nonlinear_solver|nonlinear_solver_extended|nonlinear_solver_failure_paths|operating_point|impeller_aerodynamics|impeller_prewhirl|impeller_continuity|impeller_continuity_independent|impeller_conservation|impeller_system|impeller_system_loss|incidence_loss|vaneless_diffuser|diffuser_characteristic|vaneless_diffuser_friction|vaneless_diffuser_friction_verification|vaneless_diffuser_friction_analytical|vaneless_diffuser_variable_width|vaneless_diffuser_variable_width_verification|diffuser_wall_geometry|choking_independent_validation|geometry_coupled_choking|complete_stage_reference|complete_stage_friction|complete_stage_friction_verification|complete_stage_recovery_geometry|complete_stage_reduced_recovery|impeller_diffuser_integration)$'
~~~

This is a regression subset, not an experimental gate. On numerical changes,
rerun the full suite. Prefer independent formulations and refinement studies
over tests that repeat the same code path.

### 3. Run available external comparisons

~~~sh
cmake --build build-validation --parallel 8
ctest --test-dir build-validation --output-on-failure \
  -R '^(eckardt_[ab]_external_map_screen|eckardt_b_heldout_speedline_qualification|eckardt_a_choi_kang_station_matched_screen|sandia_sco2_impeller_component_screen|schiffmann_favrat_design_point_smoke)$'
~~~

Interpret these as separate evidence:

- Eckardt A/B are map comparisons with geometry and empirical-parameter
  assumptions. The 16,000 rpm B line is the structured holdout.
- Choi–Kang is a three-point approximate digitization at specified stations;
  digitization uncertainty is material.
- Sandia compares impeller-exit P3, not the vaned-diffuser outlet. PR and
  optional HEOS runs are not full-stage CO2 qualification.
- Schiffmann–Favrat is currently a geometry/EOS/solver smoke test, not an
  accuracy pass.
- RADIVER requires authorized local data and must respect its license.

### 4. Report residuals without conflating causes

For every comparison, record source, station, measured quantity, units, inlet
state, speed/flow, geometry source, EOS, loss/slip assumptions, digitization
uncertainty, and any calibration split. Report pointwise errors, bias, MAE,
maximum absolute error, and failures; use MAPE only where meaningful. Do not
silently drop failed or out-of-domain points.

Calibration is not independent validation. Eckardt B's fit uses 10,000,
12,000, and 14,000 rpm and scores the untouched 16,000 rpm line; its
provisional 7% maximum-error gate applies only to this same-family calibrated
holdout.

### 5. Close the largest gaps

1. Resolve the disabled passage_turning_sensitivity test and rerun all tests.
2. Obtain a geometry-complete independent full-stage dataset with known
   measurement planes, inlet conditions, and uncertainty. NASA HECC access
   previously returned HTTP 403; licensed RADIVER can support diffuser work.
3. Use a reference-quality EOS in CO2 compressor calculations before assigning
   Sandia residuals to aerodynamic closures.
4. Digitize TUM vane/surge curves with uncertainty for cross-machine interface
   checks; obtain the selected unit's own maps for physical qualification.
5. Refine the multi-step recycle transient in time and compare against measured
   pipe/compressor behavior. Quantify the isothermal adapter's limits.

## Reporting rule

State exactly what a result supports: numerical identity, matched component
quantity, calibrated same-family holdout, or coupled numerical transient.
Several limited screens do not combine into broad compressor qualification.
