# Eckardt A efficiency slope screen across speed lines

Digitized efficiency points: 21. The plot symbols are approximate (about ±0.01 absolute efficiency).
Ablation rows are checked against the current digitization before analysis.
Slopes are unweighted least-squares fits of efficiency versus mass flow.
The full 14,000 rpm line is the primary result; a second row omits its 5.69 kg/s symbol only to show that point's influence.

| Scenario | Scope | rpm | n | Measured slope | Model slope | Slope error | Efficiency bias | MAE | Max abs | Within ±1 pp |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | all_digitized_points | 10000 | 6 | -1.218 | +0.104 | +1.322 | -4.26 pp | 4.26 pp | 6.44 pp | 0/6 |
| all_losses_enthalpy_equivalent | all_digitized_points | 12000 | 5 | -1.338 | +0.762 | +2.099 | -3.76 pp | 4.44 pp | 5.80 pp | 0/5 |
| all_losses_enthalpy_equivalent | all_digitized_points | 14000 | 5 | -2.466 | +0.928 | +3.394 | -2.83 pp | 5.93 pp | 7.75 pp | 0/5 |
| all_losses_enthalpy_equivalent | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.243 | -1.944 | -5.48 pp | 5.48 pp | 7.02 pp | 0/4 |
| all_losses_enthalpy_equivalent | all_digitized_points | 16000 | 5 | +2.049 | +1.088 | -0.961 | -4.17 pp | 4.17 pp | 4.93 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.536 | +0.682 | +3.58 pp | 3.58 pp | 5.99 pp | 0/6 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.325 | +1.663 | +4.49 pp | 4.49 pp | 9.39 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.615 | +3.081 | +5.72 pp | 5.72 pp | 15.88 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.976 | -2.211 | +3.18 pp | 3.18 pp | 5.52 pp | 0/4 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.883 | -1.165 | +4.73 pp | 4.73 pp | 6.51 pp | 0/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.250 | +0.968 | +0.23 pp | 1.44 pp | 2.84 pp | 2/6 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.522 | +1.860 | +0.96 pp | 1.90 pp | 6.11 pp | 2/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.758 | +3.224 | +2.06 pp | 3.59 pp | 12.41 pp | 1/5 |
| mixing_chi_0.10_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.099 | -2.089 | -0.53 pp | 1.39 pp | 2.22 pp | 1/4 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.978 | -1.071 | +0.91 pp | 0.91 pp | 2.60 pp | 3/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | +0.498 | +1.715 | -9.81 pp | 9.81 pp | 12.21 pp | 0/6 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +1.020 | +2.358 | -9.57 pp | 9.57 pp | 11.66 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +1.108 | +3.574 | -8.84 pp | 9.64 pp | 12.95 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.390 | -1.798 | -11.55 pp | 11.55 pp | 12.95 pp | 0/4 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +1.197 | -0.851 | -10.39 pp | 10.39 pp | 11.10 pp | 0/5 |
| without_skin_friction | all_digitized_points | 10000 | 6 | -1.218 | +1.163 | +2.380 | -1.49 pp | 2.92 pp | 4.20 pp | 0/6 |
| without_skin_friction | all_digitized_points | 12000 | 5 | -1.338 | +1.459 | +2.797 | -1.54 pp | 3.46 pp | 4.79 pp | 0/5 |
| without_skin_friction | all_digitized_points | 14000 | 5 | -2.466 | +1.423 | +3.889 | -0.93 pp | 5.05 pp | 10.29 pp | 0/5 |
| without_skin_friction | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.689 | -1.498 | -3.74 pp | 3.74 pp | 4.87 pp | 0/4 |
| without_skin_friction | all_digitized_points | 16000 | 5 | +2.049 | +1.430 | -0.618 | -2.60 pp | 2.60 pp | 3.32 pp | 0/5 |

## Adjacent-flow slopes for the baseline model

These secant slopes compare each neighboring pair of digitized flows. They show where the residual changes along a speed line and avoid compressing a local kink into one whole-line regression.

| rpm | Flow interval (kg/s) | Measured slope | Model slope | Slope error |
|---:|---:|---:|---:|---:|
| 10000 | 2.340–2.610 | +6.667 | +2.707 | -3.959 |
| 10000 | 2.610–3.110 | +4.000 | +1.826 | -2.174 |
| 10000 | 3.110–3.710 | -0.500 | +0.610 | +1.110 |
| 10000 | 3.710–4.280 | -8.421 | -0.759 | +7.662 |
| 10000 | 4.280–4.910 | -0.952 | -2.376 | -1.424 |
| 12000 | 2.860–3.420 | +5.179 | +2.039 | -3.140 |
| 12000 | 3.420–3.990 | +1.404 | +1.171 | -0.232 |
| 12000 | 3.990–4.540 | -0.909 | +0.354 | +1.263 |
| 12000 | 4.540–5.110 | -12.456 | -0.518 | +11.938 |
| 14000 | 3.300–3.990 | +4.928 | +1.834 | -3.094 |
| 14000 | 3.990–4.560 | +2.105 | +1.154 | -0.951 |
| 14000 | 4.560–5.110 | +2.364 | +0.570 | -1.794 |
| 14000 | 5.110–5.690 | -25.517 | -0.042 | +25.475 |
| 16000 | 3.760–4.310 | +4.182 | +1.772 | -2.410 |
| 16000 | 4.310–4.880 | +2.982 | +1.307 | -1.675 |
| 16000 | 4.880–5.440 | +0.357 | +0.863 | +0.505 |
| 16000 | 5.440–6.000 | +1.071 | +0.421 | -0.650 |

The skin-friction and chi rows are one-factor sensitivity cases; they do not establish the physical correctness of either setting.
Agreement on approximate digitized points is a model screen, not independent validation or a calibrated loss law.

Complete summary: `speedline_slope_summary.csv`; point residuals: `pointwise_efficiency_residuals.csv`; adjacent-pair slopes: `adjacent_flow_slopes.csv`.
