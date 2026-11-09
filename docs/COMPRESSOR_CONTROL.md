# Compressor operator control modes

The controller API is in \`include/turbocompgpu/compressor_control.hpp\`.
It models the supervisory process controller separately from the compressor
thermodynamics. The operator selects one primary mode and setpoint; the
controller produces a driver-speed command. Local anti-surge logic may issue
a recycle-valve command at the same time.

## Primary modes

| Mode | Operator setpoint | Feedback measurement | Manipulated variable |
|---|---|---|---|
| Discharge pressure | Discharge pressure, Pa | Discharge pressure, Pa | Driver speed, rpm |
| Suction pressure | Suction pressure, Pa | Suction pressure, Pa | Driver speed, rpm |
| Pipeline mass flow | Net flow delivered to the pipeline, kg/s | Net exported flow, kg/s | Driver speed, rpm |
| Manual speed | Driver speed, rpm | None | Driver speed, rpm |

Each automatic mode has independent PI gains and an explicit control direction.
The gains must be tuned for the compressor, driver, and connected pipeline. The
API applies minimum/maximum speed limits and a speed ramp limit. It uses
conditional integration at those limits to reduce windup. This is a simulation
control block, not a factory tuning or safety instrumented function.

## Anti-surge action

When enabled, the controller requires a surge-margin input from a separately
qualified compressor map. It drives the recycle-valve opening from 0 to 1 and
rate-limits that command. The override flag becomes active at or below the
configured target margin.

The compressor model does **not** predict a surge line or provide selected-
machine guide-vane, recycle-valve, driver-limit, or aftercooler calibration.
The component and two-pipe adapters do support guide-vane calibration,
recycle-flow accounting, and optional cooler-return callbacks when those are
supplied. Do not enable anti-surge without a margin from a validated map.
A missing margin causes the controller call to fail instead of silently
claiming surge protection. The pipe/recycle network model must apply the valve
command and conserve mass and energy.

With recycle active, compressor internal through-flow exceeds net pipeline
export flow. The controller's flow setpoint uses the latter. The compressor
physics solver must receive internal through-flow, while the connected network
accounts for the recycle branch.

## Pipe-coupling interface

The optional `transient_pipe_coupling.hpp` adapter connects compressor ports
to two transient pipe segments and supports recycle. The adapter accepts a
vane angle only when its `CompressorComponent` has a bounded, calibrated
`GuideVaneMap`. A `GuideVaneMap` maps angle, speed, flow, and inlet state to
prewhirl and entropy generation. These callbacks are interfaces; the selected
machine's calibration data are still missing.

The current pipe integration evaluates the compressor quasi-steadily at each
pipe step and holds segment temperature and compressibility fixed. It does
not apply controller PI state, shaft inertia, actuator lag, or a thermal
energy balance inside the pipe. See
`docs/TWO_PIPE_TRANSIENT_COUPLING.md` for the tested cases and commands.
