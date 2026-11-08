# Eckardt A diffuser-traverse comparison

Compared the smooth-wall, Galvas-incidence mean-line predictions with the six published 14,000 rpm mass-averaged traverse points.
Model points are evaluated using the published orifice mass-flow column. Near-station flow angle is compared only at r/r2=1.075; the paper's Table 4 does not report a far-station angle.

| Quantity | Points | Bias | Mean absolute error | Maximum absolute error |
|---|---:|---:|---:|---:|
| Near total pressure, r/r2=1.075 | 6 | -17.98% | 17.98% | 27.22% |
| Far total pressure, r/r2=1.6875 | 6 | -14.70% | 14.70% | 24.75% |
| Far static pressure, r/r2=1.6875 | 6 | -14.28% | 14.28% | 24.43% |
| Near flow angle | 6 | -7.45 deg | 7.45 deg | 17.76 deg |
| Diffuser total-pressure loss (near to far) | 6 | measured mean 4.77%; predicted mean 0.95% | mean absolute difference 3.82 percentage points | max absolute difference 4.93 percentage points |

The comparison uses integrated mass averages, not the radial profiles. The paper says this one-speed-line dataset was supplied directly by Eckardt; no equivalent traverse set is available here for other speeds. Treat this as an independent component-level check, not a general compressor qualification.

Pointwise results: `eckardt_diffuser_traverse_comparison_orifice_flow.csv`.
