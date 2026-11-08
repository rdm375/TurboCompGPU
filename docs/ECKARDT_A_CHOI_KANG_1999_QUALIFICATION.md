# Eckardt A: Choi–Kang 1999 Figure 2 screen

## Purpose and evidence level

This is a same-family comparison for backswept Eckardt Impeller A at
14,000 rpm. It uses three experimental points digitized from Figure 2 in
Choi and Kang (1999). Their paper reports the total-pressure measurement at
`r/r2=1.075`, the total-temperature measurement at `r/r2=1.687`, and an
efficiency formed using those mixed radial stations. The benchmark predicts
total pressure and total temperature at the same two radii.

The digitization is approximate. A provisional uncertainty of ±0.015 is
assigned to pressure ratio, the normalized temperature-rise coefficient, and
efficiency based on the plot resolution; the source does not publish these
as tabulated values or provide uncertainty bars for these particular points.
This is a one-speed, three-point screen, not qualification of the map or a
general accuracy claim.

## Reproduction

Build the benchmark, then run:

```sh
cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

python3 scripts/qualify_eckardt_a_choi_kang_fig2.py \
  --executable ./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep
```

The input values are preserved in
[`data/validation/eckardt_impeller_a_choi_kang_1999_fig2.csv`](../data/validation/eckardt_impeller_a_choi_kang_1999_fig2.csv).
The runner asks the C++ benchmark for predictions at the exact source flows,
checks that it used the backswept Impeller A geometry, and writes pointwise
results plus a short report under
`benchmarks/results-eckardt-impeller-a/choi-kang-1999-fig2/`.

## Result

| Flow / design flow | Flow (kg/s) | Pressure-ratio error | Temperature coefficient error | Efficiency error |
|---:|---:|---:|---:|---:|
| 0.75 | 3.405 | −2.86% | +0.0665 | −10.60 pp |
| 1.00 | 4.540 | −3.51% | +0.0542 | −10.91 pp |
| 1.50 | 6.810 | +0.92% | +0.0385 | −1.44 pp |
| **Mean absolute error** | **3 points** | **2.43%** | **0.0531** | **7.65 pp** |

Here the temperature-rise coefficient is `cp*(T03-T01)/U2^2`, using total
temperature at `r/r2=1.687`. The efficiency is recomputed from the predicted
near-plane total-pressure ratio and far-plane total temperature, matching the
source's mixed-plane convention. The model is close in near-plane pressure
ratio at the high-flow point, while it predicts a larger temperature rise at
all three points. This suggests that the design-speed energy/efficiency
closure deserves attention; it does not identify a single physical loss as
the cause. That interpretation is limited by the manual digitization and
different averaging of the original pressure and temperature traverses.

## Source

Y. S. Choi and S. H. Kang, “Application of Through-Flow Calculation to Design
and Performance Prediction of Centrifugal Compressor,” *International Journal
of Rotating Machinery*, vol. 5, no. 1, pp. 17–33 (1999), Figure 2. DOI:
[10.1155/S1023621X99000020](https://doi.org/10.1155/S1023621X99000020).
The paper describes the case as the backswept Eckardt Impeller A and states
that the plotted measured performance was obtained by spatial averaging
Eckardt's measurements. It also explains the two radial measurement stations
used for the efficiency calculation.

## Next qualification work

1. Recover tabulated values or digitize additional Figure 2 points with a
   documented pixel-to-axis procedure and uncertainty estimate.
2. Obtain and digitize the A-specific internal velocity and shroud-pressure
   profiles in the same paper; compare only variables the current model
   actually predicts at matching stations.
3. Compare the impeller specific-work/entropy closure against the measured
   temperature rise before adjusting individual loss correlations.
4. Keep the 10–16 krpm map screening separate from this one-speed component
   screen, and define acceptance tolerances before calling either qualified.
