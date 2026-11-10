# Eckardt A efficiency slope screen across speed lines

Digitized efficiency points: 21. The plot symbols are approximate (about ±0.01 absolute efficiency).
Ablation rows are checked against the current digitization before analysis.
Slopes are unweighted least-squares fits of efficiency versus mass flow.
The full 14,000 rpm line is the primary result; a second row omits its 5.69 kg/s symbol only to show that point's influence.

| Scenario | Scope | rpm | n | Measured slope | Model slope | Slope error | Efficiency bias | MAE | Max abs | Within ±1 pp |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | all_digitized_points | 10000 | 6 | -1.218 | -0.997 | +0.221 | -7.40 pp | 7.40 pp | 9.00 pp | 0/6 |
| all_losses_enthalpy_equivalent | all_digitized_points | 12000 | 5 | -1.338 | +0.078 | +1.416 | -6.35 pp | 6.35 pp | 8.29 pp | 0/5 |
| all_losses_enthalpy_equivalent | all_digitized_points | 14000 | 5 | -2.466 | +0.458 | +2.924 | -5.11 pp | 7.05 pp | 9.54 pp | 0/5 |
| all_losses_enthalpy_equivalent | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.834 | -2.353 | -7.61 pp | 7.61 pp | 9.54 pp | 0/4 |
| all_losses_enthalpy_equivalent | all_digitized_points | 16000 | 5 | +2.049 | +0.773 | -1.276 | -6.13 pp | 6.13 pp | 7.29 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -1.729 | -0.512 | +0.13 pp | 1.27 pp | 1.87 pp | 2/6 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | -0.424 | +0.914 | +1.62 pp | 1.84 pp | 5.56 pp | 3/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.096 | +2.562 | +3.18 pp | 3.75 pp | 12.66 pp | 1/5 |
| mixing_chi_0.05_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.524 | -2.664 | +0.81 pp | 1.52 pp | 3.55 pp | 1/4 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.532 | -1.517 | +2.54 pp | 2.54 pp | 4.68 pp | 0/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -1.404 | -0.187 | -3.09 pp | 3.09 pp | 4.57 pp | 0/6 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | -0.199 | +1.139 | -1.79 pp | 2.76 pp | 3.70 pp | 1/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.259 | +2.725 | -0.37 pp | 4.10 pp | 9.32 pp | 1/5 |
| mixing_chi_0.10_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.665 | -2.522 | -2.80 pp | 2.80 pp | 4.90 pp | 1/4 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.642 | -1.406 | -1.18 pp | 1.52 pp | 2.51 pp | 2/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.536 | +0.681 | -12.74 pp | 12.74 pp | 14.58 pp | 0/6 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.383 | +1.721 | -11.97 pp | 11.97 pp | 13.97 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.673 | +3.139 | -10.95 pp | 10.95 pp | 15.27 pp | 1/5 |
| mixing_chi_0.20_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.012 | -2.176 | -13.51 pp | 13.51 pp | 15.27 pp | 0/4 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.908 | -1.141 | -12.19 pp | 12.19 pp | 13.17 pp | 0/5 |
| without_skin_friction | all_digitized_points | 10000 | 6 | -1.218 | +1.163 | +2.380 | -1.49 pp | 2.92 pp | 4.20 pp | 0/6 |
| without_skin_friction | all_digitized_points | 12000 | 5 | -1.338 | +1.459 | +2.797 | -1.54 pp | 3.46 pp | 4.79 pp | 0/5 |
| without_skin_friction | all_digitized_points | 14000 | 5 | -2.466 | +1.423 | +3.889 | -0.93 pp | 5.05 pp | 10.29 pp | 0/5 |
| without_skin_friction | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.689 | -1.498 | -3.74 pp | 3.74 pp | 4.87 pp | 0/4 |
| without_skin_friction | all_digitized_points | 16000 | 5 | +2.049 | +1.430 | -0.618 | -2.60 pp | 2.60 pp | 3.32 pp | 0/5 |

## Adjacent-flow slopes for the baseline model

These secant slopes compare each neighboring pair of digitized flows. They show where the residual changes along a speed line and avoid compressing a local kink into one whole-line regression.

| rpm | Flow interval (kg/s) | Measured slope | Model slope | Slope error |
|---:|---:|---:|---:|---:|
| 10000 | 2.340–2.610 | +6.667 | +2.182 | -4.485 |
| 10000 | 2.610–3.110 | +4.000 | +1.169 | -2.831 |
| 10000 | 3.110–3.710 | -0.500 | -0.296 | +0.204 |
| 10000 | 3.710–4.280 | -8.421 | -2.056 | +6.365 |
| 10000 | 4.280–4.910 | -0.952 | -4.227 | -3.275 |
| 12000 | 2.860–3.420 | +5.179 | +1.579 | -3.600 |
| 12000 | 3.420–3.990 | +1.404 | +0.601 | -0.802 |
| 12000 | 3.990–4.540 | -0.909 | -0.390 | +0.519 |
| 12000 | 4.540–5.110 | -12.456 | -1.504 | +10.952 |
| 14000 | 3.300–3.990 | +4.928 | +1.512 | -3.416 |
| 14000 | 3.990–4.560 | +2.105 | +0.741 | -1.365 |
| 14000 | 4.560–5.110 | +2.364 | +0.049 | -2.315 |
| 14000 | 5.110–5.690 | -25.517 | -0.711 | +24.806 |
| 16000 | 3.760–4.310 | +4.182 | +1.543 | -2.639 |
| 16000 | 4.310–4.880 | +2.982 | +1.029 | -1.953 |
| 16000 | 4.880–5.440 | +0.357 | +0.523 | +0.166 |
| 16000 | 5.440–6.000 | +1.071 | +0.003 | -1.069 |

The skin-friction and chi rows are one-factor sensitivity cases; they do not establish the physical correctness of either setting.
Agreement on approximate digitized points is a model screen, not independent validation or a calibrated loss law.

Complete summary: `speedline_slope_summary.csv`; point residuals: `pointwise_efficiency_residuals.csv`; adjacent-pair slopes: `adjacent_flow_slopes.csv`.
