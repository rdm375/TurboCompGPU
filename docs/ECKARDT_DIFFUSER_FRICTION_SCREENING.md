# Eckardt A diffuser-friction screening

## Purpose

This experiment asks whether wall friction in a vaneless diffuser can explain
part of the positive pressure-ratio residual in the Eckardt Impeller A map
comparison. It is a sensitivity study, not a calibration or qualification.

## Model used

The benchmark reuses TurboCompGPU's variable-width, frictional vaneless-
diffuser integrator. The impeller result is held fixed; the diffuser conserves
stagnation enthalpy and integrates wall-friction entropy production and angular
momentum change. The reported diffuser total-pressure ratio is `p03/p02`.

The geometry uses `r3/r2=2.69`, `b3/b2=0.51`, and a linear width taper from
the effective impeller exit width to the tabulated diffuser outlet width. The
inlet effective width is 0.0255531 m, including the impeller exit blockage;
the outlet width is 0.01326 m. The linear taper and use of effective width as
the diffuser inlet width are explicit simplifications because the source does
not supply a radial width profile or station-2 blockage measurement.

The factor sweep holds a Darcy friction factor constant over each radial
integration. A second run estimates a smooth-wall Darcy factor at each point
from the diffuser-inlet Reynolds number, using Sutherland air viscosity and a
Petukhov smooth-pipe relation. That estimate is also held constant through
the diffuser. Neither approximation models jet-wake mixing, secondary flow,
separation, or a collector/volute.

## Results

All 17 points completed for each case. The no-friction row reproduces the
corrected-geometry lossless baseline. The smooth-wall estimate averages
`f_D=0.01193` over inlet Reynolds numbers averaging `8.94e5`.

| Diffuser friction case | Mean diffuser total-pressure loss | Bias | MAPE | Maximum absolute error |
|---|---:|---:|---:|---:|
| `f_D=0` | ~0% | +4.42% | 4.42% | 11.37% |
| `f_D=0.010` | 1.67% | +2.66% | 2.66% | 9.06% |
| Smooth-wall inlet-Re estimate | 1.93% | +2.39% | 2.39% | 8.78% |
| `f_D=0.020` | 3.07% | +1.18% | 1.43% | 7.04% |
| `f_D=0.030` | 4.29% | −0.10% | 1.16% | 5.26% |
| `f_D=0.040` | 5.36% | −1.23% | 1.67% | 3.66% |

Increasing wall friction reduces the overprediction, as expected. The smooth-
wall estimate improves MAPE from 4.42% to 2.39%, but leaves an 8.78% maximum
error. `f_D=0.030` gives the lowest MAPE in this sweep, but selecting it because
it best matches these same 17 points would be fitting, not independent
validation. It is not adopted as a qualified friction factor.

The result supports a limited conclusion: diffuser wall friction can explain
part of the pressure-ratio gap in this model. It does not establish that the
real diffuser has the assumed width profile or friction factor, and it does
not explain the remaining map error. The tested model also omits diffuser
mixing and separation losses, impeller losses, and collector/volute effects.

## Reproduction

Build and run the screening benchmark:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/ThermoGPU/install

cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  /tmp/eckardt-diffuser-friction.csv
```

The CSV contains per-point predicted stage, impeller, and diffuser pressure
ratios, diffuser total-pressure loss, friction factor, Reynolds number, and
relative error for each fixed-factor and smooth-wall case.
