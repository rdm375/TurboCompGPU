# Eckardt A incidence and throat-capacity screening

## Implementation

This extends the impeller-loss screen with a Galvas inlet-incidence estimate
and a local meanline throat-capacity check. It does not claim to locate the
actual three-dimensional throat or establish a compressor choke boundary.

### Inlet incidence

The correlation summarized by Meroni et al. uses

`β1,opt = atan[(A1/Ath,1) tan(β1b)]`

and

`Δh_inc = 0.5 W1² sin²(β1,opt − β1)`.

The source inlet blade angles are referenced to the meridional axis. The
benchmark converts the relative inlet velocity triangle to that same angle
convention before calculating incidence. It uses the effective annulus area
for `A1` and an idealized meanline throat estimate for `Ath,1`.

The local throat estimate uses 20 full blades, mean inlet radius, annulus span,
the mean of the published hub and shroud inlet blade angles (51.5 degrees from
the meridional direction), 2.11 mm inlet blade thickness, and zero additional
blockage. The throat-area helper expects angle from the tangential direction,
so the area calculation uses its complement (38.5 degrees):
`Z h [(2πr/Z) sin(90°−β1b) − t]`. This gives `Ath,1 ≈ 0.0303 m²`, compared
with `A1 ≈ 0.0503 m²`, and `β1,opt ≈ 58.3 degrees`. The previous implementation
used 51.5 degrees directly in the throat helper and estimated `Ath,1 ≈
0.0390 m²`; that was a convention error. The estimate remains deliberately
visible in the pointwise CSV; the actual 3D minimum throat and spanwise camber
are not available in this geometry model.

The script reports both the Galvas-incidence case and a smooth-wall case with
that incidence loss disabled. The disabled case retains the computed
incidence diagnostics but does not add the term to entropy generation.

### Throat-capacity screen

The benchmark calls the existing geometry-coupled rotating-throat solver with
the same meanline throat estimate. It assumes steady, adiabatic, isentropic,
blade-aligned flow through a locally one-dimensional throat, with normal
relative velocity sonic and zero throat-local tangential relative velocity.
It records critical mass flow, operating-flow/capacity ratio, margin, and a
flag when a sample exceeds that estimated capacity.

This is a screening calculation: the meanline area is prescribed from sparse
geometry and the model does not locate the minimum-area section. Eckardt's
published internal-flow measurements include five laser-velocimeter areas
from the inducer to the impeller discharge, but the paper identifies them as
Impeller O. Those measurements cannot independently qualify this Impeller A
throat construction. The later published Impeller A traverse figures provide
flow-field evidence, but no machine-readable throat-area or blade-to-blade
velocity table was found in the sources reviewed. A passing flag means only
that the sample is below this estimated capacity.

## Results

The run completed all 17 pressure-ratio points and all 21 digitized-efficiency
points. With the smooth-wall diffuser reference:

| Metric | Skin friction + blade loading | Plus Galvas incidence |
|---|---:|---:|
| Pressure-ratio MAPE, 17 points | 2.44% | 2.58% |
| Pressure-ratio maximum absolute error | 3.73% | 3.86% |
| Efficiency MAE, 21 points | 6.84 pp | 7.07 pp |
| Efficiency maximum absolute error | 13.88 pp | 13.71 pp |

The incidence term modestly improves the digitized efficiency MAE but slightly
worsens pressure-ratio error. The largest efficiency miss remains at 16,000 rpm
and 6.00 kg/s. Given independent, approximate efficiency digitization and the
simplified inlet throat, this does not justify making the incidence-on result
the sole reference.

The throat screen flags **0/17** pressure-map points and **0/21** efficiency
points. After converting the inlet angle to the throat helper's tangential
convention, the smallest estimated capacity margin is **24.9%** on the pressure
map and **34.8%** on the efficiency points. The minimum estimated capacity
margin is lower than before but all samples remain subcritical in this
idealized meanline screen. This does not qualify the actual 3D throat or rule
out local choking effects.

Fixed diffuser-friction factors remain sensitivity cases. For example,
`f_D=0.03` lowers efficiency MAE further but increases pressure-ratio error,
so it is not selected as a fit.

The pressure and efficiency CSVs include the estimated incidence loss, inlet
and optimum angles, throat area, critical flow, capacity utilization, margin,
and a zero/one choke-screen flag. They contain both incidence-on and
incidence-off smooth-wall rows.

## Reproduce

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  benchmarks/results-eckardt-impeller-a/impeller_incidence_choke_pressure.csv \
  data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  benchmarks/results-eckardt-impeller-a/impeller_incidence_choke_efficiency.csv
```

## Source and limits

Meroni et al., “Design of centrifugal compressors for heat pump systems,”
*Applied Energy* 232 (2018), Table 4 (Galvas incidence and impeller losses) and
Table 5 (Eckardt A geometry). [DTU repository postprint (PDF)](https://backend.orbit.dtu.dk/ws/files/154611324/Postprint_Design_of_centrifugal_compressors_for_heat_pump_systems.pdf).

The source summarizes empirical loss correlations; their use here is not
independent validation. The capacity estimate is an idealized geometric
screen. The source-reported inlet blade angles are from the meridional axis;
the throat helper's tangential-angle convention is now converted explicitly.
Eckardt's 1976 laser measurements are for the radial-discharge Impeller O, not
this backswept Impeller A. The Impeller A traverse locations and measurements
reported in Japikse (1987) are relevant for subsequent validation, but the
published figures do not supply a digitized blade-to-blade throat-area table.
Full spanwise inlet geometry, actual 3D throat location, blade camber,
clearance, and high-speed shock losses therefore remain unqualified or
unmodeled.
