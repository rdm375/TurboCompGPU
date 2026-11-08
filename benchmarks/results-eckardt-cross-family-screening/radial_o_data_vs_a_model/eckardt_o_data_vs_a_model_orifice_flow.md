# Cross-family screen: Impeller O data vs backswept Impeller A model

**Not an Impeller A qualification.** The measurements are for radial-discharge Impeller O; the executable uses backswept Impeller A geometry. These residuals cannot be interpreted as Impeller A model errors.

Compared the model with six published 14,000 rpm mass-averaged points for exploratory cross-family context only.
Model points are evaluated using the published orifice mass-flow column. Near-station flow angle is compared only at r/r2=1.075; the paper's Table 4 does not report a far-station angle.

| Quantity | Points | Bias | Mean absolute error | Maximum absolute error |
|---|---:|---:|---:|---:|
| Impeller-exit static pressure P2 | 6 | -3.76% | 3.81% | 10.60% |
| Near total pressure, r/r2=1.075 | 6 | -17.98% | 17.98% | 27.22% |
| Far total pressure, r/r2=1.6875 | 6 | -14.70% | 14.70% | 24.75% |
| Far static pressure, r/r2=1.6875 | 6 | -14.28% | 14.28% | 24.43% |
| Near flow angle | 6 | -7.45 deg | 7.45 deg | 17.76 deg |
| Diffuser total-pressure loss (near to far) | 6 | measured mean 4.77%; predicted mean 0.95% | mean absolute difference 3.82 percentage points | max absolute difference 4.93 percentage points |

The source identifies this single-speed-line dataset as radial-discharge Impeller O and says comparable public traverse data for the backswept impeller was unavailable. Do not use this cross-family screen to tune or qualify Impeller A.

Pointwise results: `eckardt_o_data_vs_a_model_orifice_flow.csv`.
