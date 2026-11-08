# TurboCompGPU Development Roadmap

Reviewed: 2026-11-06
Code status reference: default branch main, latest commit bcb7b89 (2026-11-06)

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
| Operator control modes | **Defined at interface level; integration pending** | Pressure, net-flow, and manual-speed modes need a plant/pipeline response and tuned gains |
| Variable guide vanes | **Not modeled** | No guide-vane position input or vane-dependent operating map |
| Recycle flow | **Not modeled** | No recycle branch, valve flow law, suction mixing/cooling, or exported-flow/through-flow accounting |
| Surge map and margin | **Unavailable** | Solver convergence does not predict surge; anti-surge needs an independently qualified map |
| Pipeline coupling | **Not implemented** | No end-to-end compressor boundary connected to a transient pipe model |
| ReLU surrogate | **Preliminary baseline added; not yet qualified** | A trainer fits current C++-generated samples; its structured holdout results have not yet been recorded, and it is not connected to the pipeline |
| GPU execution | **Not implemented** | Current code and reported timings are CPU based |
| Experimental validation | **Pending** | NASA HECC archive access failed with HTTP 403 in the 2026-10-13 handoff |

For supporting details, see [Project Status](PROJECT_STATUS.md) and the linked milestone reports.

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

### A2 — Add guide-vane actuation to the physical model

Where the selected machine has variable inlet guide vanes (IGVs), add vane
position as an explicit operating input. Model vane-induced prewhirl and the
associated inlet flow/incidence and loss changes using available design or
manufacturer data. Do not treat an arbitrary prewhirl value as an IGV model
without establishing that relationship.

Generate and validate operating maps across speed and vane position. Retain a
fixed-vane configuration for machines without IGVs.

**Acceptance:** vane end stops and invalid positions are handled; the fixed-vane
case recovers the existing model; representative vane/speed points agree with
independent machine data or are clearly labeled as model-only qualification.

### A3 — Model recycle flow and its thermodynamics

Add a discharge-to-suction recycle branch with a valve-opening-to-flow
relationship, line pressure losses, and any cooler used by the station.
At suction, mix fresh pipeline gas with cooled recycle using mass and energy
balances. Feed the mixed suction state and total compressor through-flow to
the compressor solver. Return net pipeline export flow separately.

**Acceptance:** zero valve opening recovers the non-recycle case; recycle mass
and energy balances close; flow tends to zero when the valve is closed; valve
limits and pressure-ratio limits are explicit; a pressure/flow controller uses
net export flow while compressor maps use internal through-flow.

### A4 — Obtain and qualify a surge-limit map

The current solver does not predict rotating stall or surge. Obtain a
speed-dependent surge line from manufacturer performance data, commissioning
tests, or other independently validated source. Include IGV position and gas
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

These items improve the reference model but do not block the first software integration prototype:

- Obtain and process the NASA HECC geometry and measurement archive; the 2026-10-13 handoff records an HTTP 403 and confirms those files were not downloaded.
- Reconstruct actual hub/shroud geometry and compare diffuser predictions with experimental data.
- Replace exploratory or prescribed loss closures with calibrated or independently supported models.
- Establish operating validity limits. The current work does not provide a predictive surge boundary; distinguish solver convergence from physical admissibility.
- Resolve the intentionally disabled test and record a fresh full regression result for the latest source tree.

## Acceleration

Keep the validated scalar CPU model as the oracle. GPU or other accelerator work is a later implementation track, justified by measured pipeline or surrogate workload. No GPU solver is currently implemented, and past CPU timings must not be presented as GPU performance.

## Working sequence

1. A0: freeze the port, flow, sensor, and control contract.
2. A1: demonstrate steady variable-speed pressure and net-flow control.
3. B0: run the ReLU baseline and independent held-out-speed-line qualification.
4. A2: add guide-vane physics and speed/vane maps for the selected machine.
5. A3: add the recycle branch with suction mixing and exported-flow accounting.
6. A4: obtain and qualify an independent surge map; anti-surge stays disabled
   until this gate passes.
7. A5: integrate the anti-surge override and unit/pipeline operating limits.
8. A6: run transient pipe–compressor–controller comparisons with refinement.
9. B1–B3: expand, qualify, and compare the surrogate through the same
   component contract.
10. A7: add station-level multiple-unit dispatch if the use case requires it.
11. Pursue experimental validation and acceleration against identified needs.
