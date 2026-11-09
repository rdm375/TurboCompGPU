# TurboCompGPU Project Status

Reviewed: 2026-11-09
Repository: [rdm375/TurboCompGPU](https://github.com/rdm375/TurboCompGPU), branch `main`
Code baseline reviewed: `df8cc301fadbeb891a97d87e83c5585563155556`

## Project goals

TurboCompGPU is pursuing two connected deliverables:

1. A compressor component that can participate in transient pipeline
   simulations.
2. A semi-realistic ReLU surrogate that recreates the sampled-operation
   interpolation approach from the 2019 NN Interpolation hackathon.

The physics model generates samples and provides a numerical reference. Its
agreement with a surrogate is separate from agreement with experimental
compressor data.

## Current status

| Area | Status | Evidence and limitation |
|---|---|---|
| Compressor physics | Implemented; numerical tests pass on the Dell | Ideal-gas and Peng–Robinson paths, impeller and complete-stage solvers, diffuser loss options, and conservation checks. Experimental qualification remains partial. |
| A0 component contract | Implemented | Static pipe-port states, flow direction, component outputs, control measurements, and failure behavior are documented. |
| A1 steady control | Implemented and tested | Manual speed, suction pressure, discharge pressure, and net pipeline-export-flow modes are available. |
| A2 guide-vane interface | Bounded tabular map loader and uncertainty fields implemented; selected-machine map missing | Complete rectilinear grids provide prewhirl, entropy generation, and conservative interpolation bounds without extrapolation. A vane angle without a bound map remains rejected. `calibrated_map_interpolation` passed on Dell 7710 on 2026-11-09. |
| A3 recycle thermodynamics | Coupling and tabular hardware-map loaders implemented; physical maps missing | Valve flow/uncertainty and flow-dependent cooler return/uncertainty grids can be loaded. Recycle mass/energy closure and two-pipe transient integration passed on Dell. `calibrated_map_interpolation` passed on Dell 7710 on 2026-11-09. |
| Two-pipe network adapter | Implemented | `transient_pipe_coupling.hpp` couples two `pipe_sim` legs to the compressor and includes recycle and delivery-pressure helpers. The pipe model is isothermal with fixed segment temperature and compressibility. |
| Multi-step recycle transient | Implemented and targeted integration test passed | The disturbance schedule carries recycle state across changing demand and valve commands and checks pipe linepack, component mass, and suction-mixer energy closure at each step. Dell CTest `two_pipe_compressor_coupling` passed in 62.31 s; time-step refinement and controller/actuator dynamics remain open. |
| Surge map / anti-surge | Bounded interpolation with uncertainty adjustment implemented; protection not qualified | A corrected-speed/vane-angle table computes surge flow, nominal margin, and uncertainty-adjusted conservative margin with strict domain checks. `calibrated_map_interpolation` passed on Dell 7710 on 2026-11-09. No selected-machine surge values are installed; protection must remain disabled without a qualified map. |
| ReLU surrogate | Training/export and leave-one-speed-line-out evaluator implemented; qualification result pending | The evaluator reports whole-line interior interpolation and endpoint extrapolation. No current qualification artifact or approved pipeline error gate is recorded. |
| Surrogate hardware | Benchmarked, with limited benefit on the small network | The earlier Dell measurements show this 6,368-MAC model is faster on the M3000M GPU than CPU, while Hailo-10H throughput is lower. Hardware throughput does not establish generalization. |
| External experimental comparison | Partial screening only | Eckardt A/B, Choi–Kang, and Sandia component screens exist. Reported screens use approximations and model-specific calibration; they are not a broad independent validation of one compressor. |

## Validation facts and limits

- The Peng–Robinson complete-stage network is 2→64→64→32→3 with 6,368
  multiply-accumulates per evaluation. Hailo validation on the 441 training
  rows is in-sample quantization validation, not held-out speed-line evidence.
- The documented surrogate qualification must hold out a complete speed line,
  fit feature and target scaling on training lines only, and report interior
  interpolation separately from endpoint extrapolation. Run it on the Dell
  after generating the current boundary grid; this checkout has no saved run
  artifact, and this environment lacks the CMake/PyTorch setup to produce one.
- Eckardt B's 16,000 rpm line was held out after fitting an empirical
  effective head/slip parameter on 10,000–14,000 rpm lines. Dell results were
  3.785% held-out MAPE and 5.151% maximum absolute error under a provisional
  7% max-error gate. This qualifies one same-family, calibrated line only.
- The Sandia sCO2 component screen with CoolProp HEOS converged on all 12
  operating points after the inlet velocity bracket fix. Raw outlet-pressure
  error was 8.65% MAPE and 16.09% maximum. HEOS spot checks agree with the
  cited Span–Wagner values; Peng–Robinson has material near-critical property
  errors. Neither is a full-stage CO2 qualification.
- On the Dell, the enabled regression suite was previously reported as 75
  passed, 0 failed, with one intentionally disabled test. After this patch,
  the targeted `two_pipe_compressor_coupling` CTest passed (62.31 s). The full
  regression suite has not yet been rerun with the extension.

## Remaining project-level work

1. Run and review the surrogate's independent whole-speed-line report, then
   choose and document a use-case error gate before any pipeline use.
2. The multi-step recycle disturbance test passed on the Dell. Compare its
   transient profile under time-step refinement, then add controller dynamics
   when actuator data are available.
3. A source-reported TUM IGV summary is now transcribed in
   `data/validation/tum_chen_2017_vigv_reported_results.csv` (Table 3.3 and
   Figure 3.17 text). Digitize the actual per-speed performance and surge
   curves with uncertainty and provenance. This is a different compressor:
   it can test interface behavior and trends, but cannot calibrate this
   machine's vane callback or authorize its anti-surge control.
4. Obtain the selected compressor's own guide-vane schedule, valve/cooler
   data, and surge line from a manufacturer, test campaign, or licensed map.
5. Run a measured pipeline/compressor transient and compare the isothermal
   adapter with an energy-transporting pipe model before operational claims.

## Related documentation

- [Development roadmap](ROADMAP.md)
- [Compressor controls](COMPRESSOR_CONTROL.md)
- [Two-pipe transient coupling](TWO_PIPE_TRANSIENT_COUPLING.md)
- [External compressor qualification data and limits](EXTERNAL_COMPRESSOR_QUALIFICATION.md)
- [Project gap-close data and run record](PROJECT_GAP_CLOSE_2026-11.md)
- [Calibrated map format and qualification requirements](CALIBRATED_MAPS.md)
- [A0–A3 implementation contract](A0_A3_IMPLEMENTATION.md)
- [Surrogate qualification workflow](../benchmarks/surrogate/README.md)
