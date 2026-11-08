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

The present compressor model does **not** predict a surge line, recycle flow,
guide-vane effects, driver power limits, or aftercooler response. Do not enable
anti-surge without supplying a surge margin derived from a validated map.
A missing margin causes the controller call to fail instead of silently
claiming surge protection. The pipe/recycle network model must apply the valve
command and conserve mass and energy.

With recycle active, compressor internal through-flow exceeds net pipeline
export flow. The controller's flow setpoint uses the latter. The compressor
physics solver must receive internal through-flow, while the connected network
accounts for the recycle branch.

## Coupling loop

At each control interval, the pipeline solver supplies measured suction and
discharge pressures and net exported flow. Call \`CompressorController::step\`
with the selected mode, target, measurements, and time step. Apply the returned
speed command to the compressor operating point and the recycle command to an
explicit recycle branch. The compressor's outlet state and shaft power then
feed back into the pipe-network update.

The current TurboCompGPU tree has no pipe-network boundary component. This
controller therefore provides the operator/control contract and actuator
commands; a separate coupling milestone must connect those commands to the
compressor and pipe equations. Fixed-speed inlet-guide-vane control is also
not represented because the current aerodynamic model has no guide-vane
position input.
