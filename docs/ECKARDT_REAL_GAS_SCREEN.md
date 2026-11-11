# Eckardt A thermodynamic property-model screen

## Purpose

The Eckardt A map runner previously instantiated `thermomodel::IdealGas` and
computed total-to-total efficiency from the constant-\(\gamma\) temperature
relation. TurboCompGPU's core solver supports Peng–Robinson, but this
validation runner was not using it.

The runner now accepts `ideal-gas` or `pr-air` as its final optional argument.
Both choices use the same measured points, geometry, loss settings, and
property-consistent thermodynamic calculations. The default remains
`ideal-gas`; its results reproduce the previous runner to floating-point
precision.

## Peng–Robinson pseudo-air assumptions

The `pr-air` option uses a Peng–Robinson pseudo-pure dry-air model with
constant ideal-gas heat capacity. The program writes a
`thermodynamic_model.txt` beside each output set with these parameters:

- Molar mass: 0.0289652 kg/mol
- Critical temperature: 132.53 K
- Critical pressure: 3.786 MPa
- Acentric factor: 0.035
- Ideal-gas heat capacity: 1004.5 J/(kg K)
- Reference state: 298.15 K, 101325 Pa

This is an EOS sensitivity screen, not a multicomponent air model; constant
heat capacity remains an approximation. Both runs retain the benchmark's
existing inlet stagnation condition of 288.15 K and 101000 Pa.

The selected property model is used throughout: inlet static-state recovery,
isentropic sound speed, rotating-throat state recovery, enthalpy-to-entropy
conversion for internal losses, external-loss enthalpy addition, diffuser
integration, and total-to-total efficiency. Efficiency is calculated as

\[
\eta_{tt} = \frac{h_{03s}-h_{01}}{h_{03}-h_{01}},
\]

where the isentropic outlet state has the inlet total entropy and predicted
outlet total pressure. With constant heat capacity and ideal-gas properties,
this is equivalent to the previous temperature-ratio expression.

## Reproduction

Build the existing target and run the cases into separate directories so the
model record for one run does not overwrite the other:

```bash
cmake --build build-model-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

exe=./build-model-validation/turbocompgpu_eckardt_diffuser_friction_sweep
input=data/validation/eckardt_impeller_a.csv
efficiency=data/validation/eckardt_impeller_a_efficiency_digitized.csv

"$exe" "$input" \
  benchmarks/results-eckardt-impeller-a/ideal-gas/pressure.csv \
  "$efficiency" \
  benchmarks/results-eckardt-impeller-a/ideal-gas/efficiency.csv \
  ideal-gas

"$exe" "$input" \
  benchmarks/results-eckardt-impeller-a/pr-air/pressure.csv \
  "$efficiency" \
  benchmarks/results-eckardt-impeller-a/pr-air/efficiency.csv \
  pr-air
```

## Current screen result

A local CPU build was syntax-checked, then run against the full efficiency
map. The pre-existing ideal-gas outputs were preserved: the largest absolute
change was about 1e-11 in pressure ratio and 1e-12 in efficiency. The
Peng–Robinson run completed all 21 efficiency points and the full pressure
map.

For the smooth-wall, Galvas-incidence case at 14,000 rpm, the PR-minus-ideal
predicted efficiency changes at 3.30, 3.99, 4.56, 5.11, and 5.69 kg/s were
-0.00439, -0.00352, -0.00240, -0.00088, and +0.00135 percentage points. The
mean absolute change was 0.00251 percentage points. The fitted efficiency
slope changed from +0.94012 to +0.94250 percentage points per kg/s.

Those changes are far smaller than the efficiency residuals and the slope
discrepancy under investigation. The Peng–Robinson correction does not explain
the observed mismatch at this low-pressure test condition. The plotted
measurements remain approximate digitizations, and the 5.69 kg/s efficiency
symbol remains a questionable high-flow point.

This screen does not close the performance-validation gap. It shows that the
ideal-gas property path in this benchmark is not the source of the larger
slope mismatch. Continue the investigation with the measurement-plane and
efficiency definition, reconstructed versus reported geometry, and
aerodynamic losses that vary with flow.
