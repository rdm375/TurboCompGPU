# Two-pipe transient compressor coupling

`include/turbocompgpu/transient_pipe_coupling.hpp` connects a compressor
component to two instances of the optional `pipe_sim::transient_step` solver.
The adapter advances the pipe states with an implicit timestep, closes
compressor through-flow against the downstream pipe's computed inlet flow, and
keeps downstream pipeline export as a separate outlet boundary.

## Recycle branch

The recycle-aware adapter advances fresh feed through the upstream pipe,
returns discharge flow through a caller-supplied valve map and return-state
model, mixes fresh and recycle streams at suction with the EOS, and advances
the downstream leg at net export flow. The compressor receives internal
through-flow; the downstream pipe carries net export. It reports recycle and
fresh flow, recycle and compressor-node mass residuals, suction mixer enthalpy
residual, and any cooler duty. Valve, return-line, and cooler models must be
calibrated for the machine and station.

## Integration test

The `two_pipe_compressor_coupling` test uses a real `CompressorComponent` and
two implicit pipe-segment steps. It checks single-step no-recycle and recycle
mass/energy closure, delivery-pressure and flow control, and carries pipe
states forward through a baseline export-demand change. The recycle transient
extension imposes a four-step demand/valve profile:

| Step | Net export (kg/s) | Recycle valve opening |
|---:|---:|---:|
| 1 | 5.0 | 0.00 |
| 2 | 4.7 | 0.25 |
| 3 | 5.3 | 0.50 |
| 4 | 5.0 | 0.00 |

Each step starts from the preceding pipe states and checks both pipe linepack
balances, compressor-node and recycle-flow closure, suction-mixer energy
closure, and `through-flow = fresh flow + recycle flow`. The test holds speed
fixed and treats valve positions as externally commanded; it does not model
controller dynamics or actuator lag. This is a software transient, not a
measured plant case. Run it on Dell before treating the new extension as
verified.

## Build and run

The pipe solver remains an optional source checkout:

```sh
cmake -S . -B build-coupled \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu" \
  -DTURBOCOMPGPU_PIPE_SIM_SOURCE_DIR="$HOME/Sync/transient-pipe-flow"
cmake --build build-coupled \
  --target turbocompgpu_two_pipe_compressor_coupling_test --parallel 8
ctest --test-dir build-coupled --output-on-failure \
  -R '^two_pipe_compressor_coupling$'
```

## Scope

- Pipe segments use fixed temperature and compressibility factor; the adapter
  does not transport compressor outlet enthalpy through a thermal pipe model.
- Each compressor evaluation is quasi-steady. The test does not add shaft
  inertia, controller sample/hold behavior, valve lag, or measured disturbances.
- A machine-specific IGV map and recycle-valve/cooler calibration are still
  required. The optional guide-vane input is rejected if a calibrated map was
  not supplied.
- No selected-machine surge line has been qualified. Anti-surge protection
  must remain unavailable without a valid, in-domain margin provider.
- The adapter and its synthetic tests establish software coupling and
  conservation only, not experimental compressor or pipeline accuracy.
