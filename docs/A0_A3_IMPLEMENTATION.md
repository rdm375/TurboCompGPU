# A0–A3 compressor/pipeline implementation contract

The header-only interfaces in \`include/turbocompgpu/pipeline_compressor.hpp\`
and \`include/turbocompgpu/steady_pipeline_control.hpp\` establish the
compressor boundary and the steady control/recycle calculations. They are
designed to be bound to a pipe-network evaluator; this repository does not
currently contain that pipe solver.

## A0 — Component ports and conventions

- Pipe ports use static temperature and absolute static pressure.
- Positive pipeline export is upstream-to-downstream.
- The compressor component receives **compressor through-flow** and returns
  station-3 static and total states plus shaft power.
- Guide-vane angle is optional. Supplying it requires a calibrated,
  range-bounded \`GuideVaneMap\`.
- The fixed-vane path uses explicit fixed prewhirl and entropy-generation
  inputs. Residual scales are numerical normalization values and should be
  set for the chosen gas and operating envelope.

## A1 — Steady setpoint coupling

`solve_steady_speed_setpoint` evaluates a caller-supplied plant callback at
speed limits, then bisects speed against suction pressure, discharge pressure,
or pipeline-export-flow setpoint. Manual-speed mode remains available. The
configured minimum and maximum speed are the feasible actuator bounds. The
callback must solve the connected pipe/compressor network and report its pipe
tap measurements. When a maximum shaft-power limit is configured, the callback
must also report shaft power; the power-limited search assumes power is
nondecreasing over the configured speed interval. Unreachable targets report
speed or power saturation.

## A2 — Guide-vane input

The existing impeller solver already accepts inlet prewhirl and entropy
generation. The new component adapter maps a vane command to those physical
inputs through a bounded callback. The callback must be populated from
machine/vendor/commissioning data; no universal vane-angle correlation is
invented here. Out-of-range queries fail explicitly.

## A3 — Recycle accounting

\`solve_steady_recycle\` enforces the bookkeeping identity

\`compressor through-flow = pipeline export flow + recycle flow\`.

It mixes fresh and returned gas at the suction pressure using mass-weighted
enthalpy and obtains the mixed state from the selected EOS. Valve flow comes
from a bounded calibration map. The caller supplies a return-state function
that represents throttling and any aftercooler; the provided
\`isenthalpic_recycle_return\` helper covers throttling without cooling.
A fixed-point iteration couples valve flow to compressor discharge pressure.

The reported cooler duty is positive when the return stream is cooled. A
nonconverged recycle loop is returned as a failure. The valve map and
aftercooler still require system-specific data. The permanent
`pipeline_compressor_control` CTest covers the component adapter, all
setpoint modes, power limiting, calibration-map bounds, EOS mixing, and
recycle-loop mass/energy balances. Its fixtures use deterministic callbacks,
not a calibrated compressor map.

The optional `two_pipe_compressor_coupling` CTest advances the native
`pipe_sim` solver on both sides of a real compressor component. Its inner
closure equates compressor through-flow with the downstream pipe's computed
inlet flow; pipeline export remains the downstream outlet-flow boundary.
The test applies discharge-pressure speed control, then advances further
steps after changing export demand. See
`docs/TWO_PIPE_TRANSIENT_COUPLING.md` for build commands and the isothermal,
closed-recycle scope.

## Remaining qualification needed

This two-pipe coupling does not constitute a measured pipeline case or a
fully coupled compressor/recycle/thermal transient model. Supply machine
guide-vane and recycle-valve data, extend the pipe state to carry the needed
thermal/EOS variables, add transient recycle coupling, and run the qualification
on the Dell 7710. The surge-map and anti-surge qualification remains a later
milestone.
