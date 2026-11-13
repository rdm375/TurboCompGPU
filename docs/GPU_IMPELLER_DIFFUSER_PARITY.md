# CUDA impeller and vaneless diffuser parity

The CUDA batch backend now carries each converged impeller point through the
same constant-width, zero-blockage, frictional vaneless diffuser equations used
by the CPU analytical integrator. The diffuser inlet is the impeller exit; its
outlet radius and radial RK4 step count are explicit batch inputs.

The impeller passage still uses a lumped inlet-Re Colebrook factor. In the
diffuser, every RK4 stage evaluates local PR density and velocity, computes
ISO 20765-5 LBC viscosity from local temperature, pressure, and PR density,
computes Re using `Dh = 2 b`, and evaluates Colebrook on that local Reynolds
number. Relative roughness is held constant. CPU/GPU parity independently
evaluates this closure on CPU and GPU. At 512 radial steps, each operating point
performs 2,048 local viscosity and Colebrook evaluations, in addition to the
repeated PR-mixture state evaluations.

The CPU/GPU parity example independently evaluates inlet Reynolds number,
Colebrook friction, and the lumped impeller passage entropy on CPU, compares
those loss quantities with the GPU values, and uses the CPU-computed passage
entropy in the CPU compressor solve. It then compares impeller outlet state,
Euler work, and torque, followed by diffuser outlet pressure, temperature,
density, radial and tangential velocity, stagnation-enthalpy residual, entropy
change and production, entropy balance, wall-torque integral, angular-momentum
balance, and mass residual. The parity driver covers a 5-by-4 natural-gas
screening grid: 10,000–18,000 rpm and 3–6 kg/s at fixed inlet conditions and
geometry. It checks each pair of CPU/GPU solutions and fails with the speed and
flow attached to any mismatch.

The current implementation requires zero impeller-exit blockage and a
constant-width diffuser. The smoke/parity inputs use an outlet radius of
0.320 m and 512 RK4 radial steps. The LBC implementation uses PR density in
place of the GERG-2008 density expected by ISO 20765-5, so it is a screening
model and not a claim of ISO compliance. It enforces 253.15–313.15 K and
0.1–10 MPa. Impeller passage friction also remains lumped until the API carries
internal impeller station geometry and states. Details are in
[GPU_FRICTION_CLOSURE.md](GPU_FRICTION_CLOSURE.md).

With the CUDA and dependency prefixes configured, build and run:

```sh
cmake --build build-cuda-gcc13 \
  --target turbocompgpu_gpu_impeller_batch \
           turbocompgpu_gpu_impeller_cpu_parity --parallel 8

./build-cuda-gcc13/turbocompgpu_gpu_impeller_batch
./build-cuda-gcc13/turbocompgpu_gpu_impeller_cpu_parity \
  > benchmarks/results-gpu-impeller-throughput/operating-map-parity.csv \
  2> benchmarks/results-gpu-impeller-throughput/operating-map-parity.log
cat benchmarks/results-gpu-impeller-throughput/operating-map-parity.log
```

A successful parity run ends with:

```text
CPU/GPU impeller, Colebrook, and frictional diffuser parity passed for 20 points across 5 speeds and 4 flow rates.
```

## Batch throughput

The throughput executable checks CPU/GPU parity over each batch before
reporting timings. It reports serial CPU, CPU-parallel, and GPU end-to-end
points per second. Device input/output buffers are allocated once when the
solver is constructed and reused; timings therefore exclude one-time setup but
include input/output copies and the synchronization required to return results.
The CPU-parallel path distributes independent operating points over OpenMP
threads when available, or a small `std::thread` worker set otherwise.

```sh
./build-cuda-gcc13/turbocompgpu_gpu_impeller_batch_throughput \
  --sizes 1,4,16,64,256 --repeats 3
```

Use `--cpu-threads N` to choose the host worker count and `--sizes` to extend or
reduce the batch sweep. Run the same command on each GPU host; compare the
reported GPU points/s and speedups, while keeping compiler, power mode, and
background load consistent. These numbers measure this screening stage and
API, not a full compressor map or a GPU-only kernel ceiling.
