# Eckardt A inlet throat and diffuser measurement alignment

## What changed

The local inlet-throat area had a blade-angle convention error. Meroni's
Eckardt A inlet angles are reported from the meridional direction, while
`ThroatGeometry` defines its blade angle from the tangential direction. The
benchmark had passed their arithmetic mean (51.5 degrees) directly to the
throat-area helper. The corrected projection uses its complement (38.5
degrees from tangential):

```text
A_th = Z h [(2 pi r / Z) sin(90 deg - beta_meridional) - t]
```

For the selected meanline inputs (`Z=20`, `r=0.1077 m`, `h=0.080 m`,
`t=2.11 mm`, and zero extra blockage), this changes `A_th` from about
`0.0390 m^2` to `0.0303 m^2`. The inlet annulus area remains about
`0.0503 m^2`. The incidence correlation continues to use 51.5 degrees in the
source's meridional convention; it does not use the converted throat-helper
angle.

## What the measurements support

The measured geometry and internal-flow evidence are not all from the same
Eckardt configuration:

- Meroni's Eckardt A geometry gives a 120 mm inlet hub diameter, 280 mm inlet
  shroud diameter, 20 full blades, 40-degree hub and 63-degree shroud inlet
  blade angles, and 2.11 mm inlet blade thickness. The benchmark uses a
  simplified meanline through these values; it does not reconstruct the
  spanwise camber or locate the minimum three-dimensional passage area.
- Eckardt's 1976 laser-velocimeter paper presents five internal-flow
  measurement areas from inducer inlet to impeller discharge, but that paper
  describes the radial-discharge Impeller O. It is useful context for the
  measurement method and flow physics, but it is not direct validation of this
  backswept Impeller A throat.
- Japikse's 1987 assessment includes the O, A, and B configurations. It reports
  total-pressure and flow-angle traverses at four diffuser radii, including
  `r/r2=1.075` nearest the impeller tip and `r/r2=1.6875` farthest outward,
  together with wall static-pressure measurements. Its discussion also
  identifies static taps at the impeller inlet, along the shroud, at the
  impeller exit, and through the vaneless diffuser. The source describes a
  large manually transcribed traverse database, but the machine-readable
  pointwise measurements were not found in the reviewed public copy.

These sources therefore provide enough evidence to correct the angle
convention and define future pointwise comparison locations. They do not
currently provide a directly comparable, digitized Eckardt A blade-to-blade
throat area or velocity profile for qualifying the 3D throat estimate.

## Diffuser station alignment

The current benchmark integrates a linearly tapered diffuser to `r3/r2=2.69`
for its Eckardt A map prediction. For `r2=0.200 m`, that model endpoint is
`r3=0.538 m`. The public detailed traverse description gives the nearest and
farthest traverse stations below; these are measurement planes, and the far
traverse is not identified as the diffuser outlet.

| Location | `r/r2` | Radius at `r2=0.200 m` | Width in current linear-taper model |
|---|---:|---:|---:|
| Near traverse | 1.075 | 0.2150 m | 0.02501 m |
| Far traverse | 1.6875 | 0.3375 m | 0.02055 m |
| Model outlet | 2.69 | 0.5380 m | 0.01326 m |

The widths at the traverse planes are model interpolation values, not
published measurements. The benchmark now writes
`diffuser_traverse_predictions.csv` with the modeled width, bulk static
pressure, total pressure, total temperature, flow angle from the meridional
plane, and both velocity components at each of these radii for every pressure
map point. The total-pressure and angle predictions are axisymmetric
meanline states; bulk static pressure is not the wall-tap pressure. The
experimental profiles should be compared at their actual traverse radii; the
far-traverse data should not be treated as an outlet pressure measurement.
Before qualifying diffuser friction, we need pointwise Impeller A traverse data
or a documented digitization of the published profiles, plus confirmation of
the pressure-map outlet plane used by the source map.

## Updated screening result

The corrected model still solves all 17 pressure-map points and 21 approximate
efficiency points. The pressure-map screening MAPE is 2.58%; the approximate
efficiency digitization MAE is 7.07 percentage points. The smallest estimated
capacity margin is 24.88% on the pressure rows and 34.81% on the efficiency
rows. The maximum calculated throat Mach is 0.457. These are model outputs
under the simplified meanline throat assumption, not independent validation
against measured throat capacity.

The angle correction changes the choke-screen result more than the overall
map-fit metrics. The new traverse CSV creates a radius-matched comparison
surface for digitized diffuser measurements, but no measurement values have
been entered yet. The throat remains unqualified: a geometric projection
correction fixes internal consistency, while the physical throat area still
needs an independent geometric or measured-flow anchor.

## Reproduction

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  benchmarks/results-eckardt-impeller-a/angle_convention_fix_pressure.csv \
  data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  benchmarks/results-eckardt-impeller-a/angle_convention_fix_efficiency.csv

python3 scripts/audit_eckardt_skin_friction.py \
  benchmarks/results-eckardt-impeller-a/loss_term_ablation.csv
```

## Sources

- Meroni, *Design and Optimization of Turbomachinery for Thermodynamic Cycles
  Utilizing Low-Temperature Heat Sources*, DTU thesis (2018), Eckardt A
  geometry in Appendix C. [DTU repository PDF](https://orbit.dtu.dk/en/publications/design-and-optimization-of-turbomachinery-for-thermodynamic-cycles-).
- Eckardt, “Detailed Flow Investigations Within a High-Speed Centrifugal
  Compressor Impeller,” *Journal of Fluids Engineering* 98(3), 390–399 (1976),
  DOI [10.1115/1.3448334](https://doi.org/10.1115/1.3448334).
- Japikse, “A Critical Evaluation of Three Centrifugal Compressors With
  Pedigree Data Sets: Part 5—Studies in Component Performance,” *Journal of
  Turbomachinery* 109(1), 1–9 (1987), DOI
  [10.1115/1.3262064](https://doi.org/10.1115/1.3262064).
