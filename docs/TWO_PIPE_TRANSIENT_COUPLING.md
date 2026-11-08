# Two-pipe transient compressor coupling

`transient_pipe_coupling.hpp` connects the compressor component to two
instances of the `pipe_sim::transient_step` solver. One implicit timestep is
advanced in this order:

1. Advance the upstream pipe with source pressure and compressor through-flow
   as its outlet flow; the resulting last-node pressure is compressor suction.
2. Evaluate the compressor at that suction state, trial through-flow, and
   selected speed.
3. Advance the downstream pipe with compressor discharge static pressure as
   its inlet pressure and pipeline export as its outlet flow. Its calculated
   inlet flow is matched to compressor through-flow by a scalar flow-closure
   solve. This is distinct from export flow while linepack changes.

Each returned segment state retains both the computed inlet flow and the
imposed outlet boundary flow. That prior boundary value is used on the next
step; it is not inferred from the last interior staggered-flow unknown.

The `advance_controlled_two_pipe_compressor` helper lets the existing A1 speed
controller evaluate that network step at trial speeds. Every trial starts from
the same old pipe states, and only the selected final step is returned for the
caller to commit.

## Build and run the integration test

The pipe solver remains an optional source checkout. Configure TurboCompGPU
with paths to ThermoGPU and the native C++ pipe solver:

```sh
cmake -S . -B build-coupled \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu" \
  -DTURBOCOMPGPU_PIPE_SIM_SOURCE_DIR="$HOME/Sync/transient-pipe-flow"
cmake --build build-coupled --target turbocompgpu_two_pipe_compressor_coupling_test --parallel 8
ctest --test-dir build-coupled --output-on-failure -R '^two_pipe_compressor_coupling$'
```

The test uses a real `CompressorComponent`, two real implicit pipe-segment
steps, and discharge-pressure speed control. It checks the compressor
discharge boundary, the compressor-to-downstream-pipe mass balance, and three
follow-on timesteps that carry the updated pipe states forward after an export
demand change. For each pipe leg it checks the theta-method linepack balance
against the old and new boundary flows. The compressor-flow search bounds
must cover the physical solution for the chosen machine and pipe states.

## Scope and limits

- The coupled operation is a single implicit timestep; callers advance their
  stored upstream and downstream states by assigning the returned states.
- The pipe solver is isothermal with fixed temperature and compressibility
  factor per segment. The adapter passes interface static pressures but does
  not transport compressor outlet enthalpy or temperature into the pipe.
- Recycle is closed in this first coupled adapter. The A3 recycle mixer and
  valve loop remain available separately, but do not yet participate in this
  transient network step.
- There is no surge-map or anti-surge control qualification.
- The present integration test is a software-coupling check, not a calibrated
  station or measured pipeline validation.
