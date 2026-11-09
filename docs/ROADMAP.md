# TurboCompGPU Development Roadmap

Reviewed: 2026-11-09
Code status reference: default branch main, latest reviewed code commit 3980363d090d3b532cf50d437f8b90cfc9a0fa21

## Project goals

The project has two intended outputs:

1. **Pipeline component:** a compressor model that can be connected to a transient pipe solver and participate in pipeline simulation.
2. **Semi-realistic surrogate candidate:** a fast ReLU interpolant of sampled compressor operation, continuing the approach explored in the original NN Interpolation hackathon.

These goals share data and interfaces, but they are not the same deliverable. The compressor physics solver is the reference model and possible sample generator. The surrogate is a fitted approximation with its own accuracy, domain, and validation requirements.

## Current status

| Area | Status | Evidence and limitation |
|---|---|---|
| Thermodynamics and nonlinear compressor solver | Implemented and numerically tested | IdealGas and PengRobinson paths, conservation checks, reference-invariance and solver validation |
| Operating points and maps | Implemented for tested domains | Deterministic sweeps and continuation; tested convergence is not proof of physical admissibility |
| Impeller/diffuser complete-stage model | Implemented with numerical verification | Frictional and variable-width diffuser work is documented; empirical calibration is still missing |
| Incidence closure | Experimental model extension | 105 paired points across three thermodynamic configurations; coefficient not experimentally calibrated |
| Diffuser profile geometry | Infrastructure present | PCHIP hub/shroud profiles and tests; actual NASA HECC hot-wall geometry is not reconstructed |
| Operator control modes | **Implemented; gains remain application-specific** | Manual speed, suction/discharge pressure, and net-export-flow modes are available. Pipe-coupled adapters and saturation behavior are tested. |
| Variable guide vanes | **Bounded map loader tested; selected-machine data missing** | `calibrated_map_interpolation` passed on Dell. Rectilinear interpolation supplies prewhirl, entropy generation, and uncertainty bounds without extrapolation. |
| Recycle flow | **Coupling and bounded hardware-map loaders tested; selected-unit data missing** | `calibrated_map_interpolation` passed on Dell. Valve flow and flow-dependent cooler return can be loaded from complete grids; hardware calibration remains absent. |
| Surge map and margin | **Bounded interpolation tested; no qualified machine map** | `calibrated_map_interpolation` passed on Dell. Anti-surge still requires a qualified selected-unit map; keep protection disabled without one. |
| Pipeline coupling | **Two-pipe adapter implemented; transient qualification incomplete** | Optional `pipe_sim` integration supports recycle; each segment currently uses fixed temperature and compressibility. |
| ReLU surrogate | **Training/export and leave-one-speed-line-out workflow implemented; qualification result pending** | Whole-line interior interpolation and endpoint extrapolation are reported separately. No saved current report or pipeline error gate is recorded. |
| GPU execution | **Not implemented** | Current code and reported timings are CPU based |
| Experimental validation | **Partial screens; machine hardware maps still missing** | Eckardt A/B, Choi–Kang, and Sandia screens exist. TUM source-reported IGV summaries are transcribed, while pointwise speed-line and surge coordinates still need digitization. TUM hardware differs from the current rotor. NASA HECC archive access failed with HTTP 403 in the 2026-10-13 handoff. |

For supporting details, see [Project Status](PROJECT_STATUS.md) and the linked milestone reports.
The current held-out-surrogate, transient-test, and physical-data acquisition
record is in [Project gap-close data and run record](PROJECT_GAP_CLOSE_2026-11.md).

The active priority is component-wise validation of the numerical compressor
model. See [Numerical model validation](NUMERICAL_MODEL_VALIDATION.md) for the
campaign and evidence boundaries. Surrogate qualification is deferred until
the numerical reference model's evidence and validity domain are reviewed.

## Workstream A — Connect a compressor to a pipeline

The operator selects a primary process objective; the control layer changes a
machine actuator to pursue it. The network determines the resulting pressure
and flow together with the compressor. A compressor operating-point solver
alone cannot close this loop.

### A0 — Freeze the unit and pipeline-port contract

Specify inlet and outlet measurement locations, pressure conventions, mass-flow
signs, thermodynamic states, control modes, actuator limits, and failure
behavior. Distinguish:

- **Pipeline export flow:** net flow delivered downstream; this is the
  measurement for a station flow setpoint.
- **Compressor through-flow:** flow through the impeller and diffuser; this is
  the mass flow supplied to the compressor physics model.
- **Recycle flow:** gas diverted from compressor discharge back to suction.
  With recycle, compressor through-flow exceeds pipeline export flow.

Use static pressure at defined pipe taps for control feedback. Keep those
sensor values distinct from total pressures returned by compressor-map
postprocessing.

**Acceptance:** a documented component/control data contract with units,
measurement locations, flow sign, thermodynamic state convention, actuator
configuration, and explicit invalid/out-of-domain behavior.

### A1 — Couple a variable-speed compressor in steady operation

Connect the compressor between upstream and downstream pipe sections. Use the
network solution to determine mass flow and boundary pressures. Implement
discharge-pressure, suction-pressure, and pipeline-export-flow control modes;
the controller adjusts shaft speed within driver and compressor limits.
Manual-speed mode remains available for fixed-command studies.

**Acceptance:** one deterministic steady pipe–compressor–pipe case for each
automatic setpoint mode; mass and energy close at the component boundary;
setpoints track when feasible; saturated commands report that the target is
unreachable within the configured speed/power envelope.

### A2 — Qualify guide-vane actuation against physical data

The component accepts a vane angle only when a bounded `GuideVaneMap` is
supplied. A rectilinear CSV loader can now bind measured prewhirl, entropy,
and conservative uncertainty bounds to the component. The remaining work is
to populate and validate that map for the selected machine. An open TUM dissertation reports stage maps at IGV settings
−20°, 0°, 20°, 40°, and 60° for three vane designs, plus surge tests. Table
3.3 and the numerical comparisons stated in the Figure 3.17 discussion are
transcribed in `data/validation/tum_chen_2017_vigv_reported_results.csv`.
Per-speed curve coordinates and surge-line digitization are still outstanding.
Its compressor differs from the project's rotor: use it to check data import,
interpolation, uncertainty, and software behavior, not to calibrate another
machine's vane map.

Generate and validate operating maps across speed and vane position. Retain a
fixed-vane configuration for machines without IGVs.

**Acceptance:** vane end stops and invalid positions are handled; the fixed-vane
case recovers the existing model; representative vane/speed points agree with
independent data for the selected compressor. Cross-machine literature screens
remain labeled as software/data-pipeline checks.

### A3 — Model recycle flow and its thermodynamics

The steady recycle solver and two-pipe adapter provide a discharge-to-suction
branch, suction mixing, and separate fresh, recycle, compressor-through-flow,
and net-export values. Complete-grid loaders now cover measured valve flow and
a flow-dependent cooler return state. Physical valve and cooler tables remain
caller-supplied until selected-unit calibration data are available. A multi-step demand/valve disturbance test
has been added after the existing single-step mass and energy checks; the Dell
integration test passed on 2026-11-09. Time-step refinement remains outstanding.

**Acceptance:** zero valve opening recovers the non-recycle case; recycle mass
and energy balances close at every transient step; flow tends to zero when the
valve is closed; valve limits are explicit; flow control uses net export while
the compressor map uses internal through-flow; time-step refinement is
reported.

### A4 — Obtain and qualify a surge-limit map

The table loader and interpolation path are implemented in
`include/turbocompgpu/calibrated_maps.hpp`; see `docs/CALIBRATED_MAPS.md` and
`data/validation/surge_line_map_template.csv`. These establish software
behavior only, not a qualified protection boundary.

The current solver does not predict rotating stall or surge. A bounded
rectilinear map can now interpolate corrected surge flow over corrected speed
and vane position and provide a margin to the controller. Obtain a
speed-dependent surge line for the selected unit from manufacturer performance
data, commissioning tests, or another independently validated source. An open
TUM dissertation provides a different-machine vane/surge dataset useful for
checking data processing and map interfaces; it cannot qualify this machine's
surge limit. Include IGV position and gas
correction where the source requires them. Define the corrected flow/head
coordinates, surge margin, uncertainty, interpolation range, and behavior
outside that range. Keep choke and driver limits distinct from the surge line.

Do not infer a surge limit from solver convergence or from the current 441-row
training grid alone.

**Acceptance:** source and units are recorded; map interpolation is checked
against withheld map points; every supported query has a bounded margin and
uncertainty; out-of-map queries return an explicit unavailable/invalid result.
Until this gate passes, anti-surge control remains unavailable.

### A5 — Integrate anti-surge override and operating limits

Use the qualified margin provider as a separate fast protection loop. It opens
the recycle valve as the operating point approaches the configured control
line, with suitable margin, rate limits, and actuator response. Anti-surge
action takes priority over the normal pressure/flow objective. Add explicit
driver speed and power limits, maximum discharge pressure, and discharge
temperature limits from the selected unit and downstream pipeline
specification. Treat shutdown trips as separate protection logic.

**Acceptance:** disturbance and setpoint-change scenarios preserve the
qualified surge margin; override priority, valve travel/rate limits, saturation,
and missing-map behavior are exercised; primary-loop windup is controlled
during override. No test may claim surge protection with an absent or
out-of-domain map.

### A6 — Validate transient pipe–compressor–controller coupling

Run the same unit with the pipe solver in transient operation. Define the
controller sample period, actuator lag/ramp response, and whether the
compressor is modeled quasi-steadily at each pipe time step. Exercise changes
in demand, pressure setpoint, speed/IGV, and recycle action. Compare time-step
refinement and coupling iteration.

**Acceptance:** repeatable transient traces report pressure/flow tracking,
speed and valve commands, surge margin when qualified, mass/energy residuals,
controller saturation, and solver failures. The steady case remains a
regression reference.

### A7 — Expand to station behavior only after the single unit

Add parallel-unit load sharing, unit starts/stops, bypasses, and station
dispatch only when the single-compressor coupling is stable. Keep station
dispatch separate from the unit-level pressure/flow and anti-surge loops.

## Workstream B — Reproduce and qualify the ReLU interpolation idea

The original NN Interpolation hackathon used ReLU networks as interpolation functions over sampled compressor-operation input/output data. It treated the compressor simulator as a source of samples rather than requiring access to the simulator internals. Preserve that as the baseline method and record the provenance and units of every sample set.

### B0 — Establish a model-generated ReLU baseline

The original hackathon artifacts are not currently available. Use the existing compressor-map samples to recreate the core ReLU interpolation method, while labeling the data as outputs of the current C++ model rather than legacy-simulator or experimental truth. The initial script is documented in RELU_SURROGATE_BASELINE.md.

**Acceptance:** run the training script, save its structured speed-line holdout metrics, and confirm that the model checkpoints can reproduce predictions from the report.

### B1 — Define a semi-realistic training domain

Use the physics-based compressor model to generate additional samples only within the model’s documented validity domain. Where legacy simulator samples exist, compare the two sources at matched conditions instead of silently treating either as experimental truth.

Vary relevant operating inputs (for example inlet state, mass flow, and speed) while preserving geometry and gas composition metadata.

**Acceptance:** a dataset manifest, domain coverage report, solver-failure counts, and a clear label for each source: legacy simulator, current physics model, or experiment.

### B2 — Validate the interpolant

Compare held-out predictions against their sample source and against independent physics-solver evaluations. Report maximum, percentile, and RMS errors per output, boundary behavior, and out-of-domain behavior. Check physical consistency, including positive absolute pressure and nonnegative power where applicable.

**Acceptance:** a qualified domain and explicit error gate chosen for the pipeline use case. Do not use random train/test splits alone when nearby operating points can leak across the split; hold out structured regions or complete speed lines as well.

### B3 — Use the surrogate in a pipeline case

Expose the surrogate through the same component contract as the reference solver. Compare coupled pipeline outputs against the reference-model coupling over the qualified domain.

**Acceptance:** repeatable pipeline-level comparisons and a defined response when queries leave the surrogate’s qualified domain.

## Model credibility work

This pass focuses on the numerical compressor model. Keep these claims
distinct:

- Numerical verification checks equations, derivatives, conservation,
  independent formulations, tolerances, and solver failure behavior.
- Component validation compares the same measured quantity at matched
  operating conditions and measurement planes.
- Calibration fits declared training data; validation scores untouched data.
- Pipe integration checks network closure and transient behavior, but does not
  establish compressor accuracy against a test rig.

The current evidence does not include a geometry-complete independent
full-stage dataset. Follow the component sequence in
[Numerical model validation](NUMERICAL_MODEL_VALIDATION.md), and do not turn
limited or cross-family screens into a broad compressor qualification.

Remaining data limits:

- NASA HECC archive access returned HTTP 403 in the 2026-10-13 handoff.
- Eckardt A/B, Choi–Kang, and Sandia use different assumptions and support
  component/performance screens only.
- RADIVER requires authorized local data and must remain subject to its signed
  license.
- Selected-machine guide-vane, recycle hardware, and surge maps are missing.
- The passage_turning_sensitivity test remains disabled; report that status
  and rerun the complete suite after resolving it.

## Acceleration

Keep the validated scalar CPU model as the oracle. GPU or other accelerator work is a later implementation track, justified by measured pipeline or surrogate workload. No GPU solver is currently implemented, and past CPU timings must not be presented as GPU performance.

## Working sequence

1. Run all enabled CTest tests and record the commit, compiler, build/EOS
   options, optional dependencies, disabled tests, and complete output.
2. Complete numerical verification by component: thermodynamic states and
   derivatives; nonlinear solver paths; impeller continuity, work, and
   conservation; diffuser friction and geometry; choking; and complete-stage
   integration.
3. Run matched external screens for Eckardt A/B, Choi–Kang, and Sandia. Use
   CoolProp HEOS when available; keep RADIVER local under its license and
   Schiffmann–Favrat labeled as a smoke test until its model inputs are
   qualified.
4. Report pointwise residuals, bias, MAE, maximum error, failures, measurement
   uncertainty, station definitions, geometry, and EOS assumptions. Do not
   combine separate datasets into one pass/fail score.
5. Keep the Eckardt B 16,000 rpm line as its documented calibrated holdout.
   Do not describe its fitted result as an uncalibrated or cross-family pass.
6. Seek a geometry-complete independent full-stage dataset; revisit NASA HECC
   access and pursue authorized RADIVER data.
7. Resolve the disabled sensitivity test, then refine the multi-step recycle
   pipe transient and compare it with measured behavior.
8. Obtain selected-machine vane, recycle, and surge data before qualifying
   hardware behavior or anti-surge protection.
9. Resume surrogate qualification after reviewing the numerical model's
   component evidence and validity domain.
10. Add station dispatch and acceleration only where the validated application
    requires them.
