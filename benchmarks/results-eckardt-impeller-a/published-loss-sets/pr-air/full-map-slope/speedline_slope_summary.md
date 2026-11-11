# Eckardt A efficiency slope screen across speed lines

Digitized efficiency points: 21. The plot symbols are approximate (about ±0.01 absolute efficiency).
Ablation rows are checked against the current digitization before analysis.
Slopes are unweighted least-squares fits of efficiency versus mass flow.
The full 14,000 rpm line is the primary result; a second row omits its 5.69 kg/s symbol only to show that point's influence.

| Scenario | Scope | rpm | n | Measured slope | Model slope | Slope error | Efficiency bias | MAE | Max abs | Within ±1 pp |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | all_digitized_points | 10000 | 6 | -1.218 | +0.143 | +1.361 | -4.35 pp | 4.35 pp | 6.55 pp | 0/6 |
| all_losses_enthalpy_equivalent | all_digitized_points | 12000 | 5 | -1.338 | +0.785 | +2.122 | -3.86 pp | 4.51 pp | 5.91 pp | 0/5 |
| all_losses_enthalpy_equivalent | all_digitized_points | 14000 | 5 | -2.466 | +0.942 | +3.409 | -2.94 pp | 6.01 pp | 7.66 pp | 0/5 |
| all_losses_enthalpy_equivalent | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.255 | -1.933 | -5.59 pp | 5.59 pp | 7.13 pp | 0/4 |
| all_losses_enthalpy_equivalent | all_digitized_points | 16000 | 5 | +2.049 | +1.095 | -0.953 | -4.28 pp | 4.28 pp | 5.04 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.492 | +0.726 | +3.49 pp | 3.49 pp | 5.92 pp | 0/6 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.352 | +1.689 | +4.37 pp | 4.37 pp | 9.31 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.631 | +3.098 | +5.60 pp | 5.60 pp | 15.79 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.989 | -2.198 | +3.05 pp | 3.05 pp | 5.38 pp | 0/4 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.892 | -1.157 | +4.60 pp | 4.60 pp | 6.38 pp | 0/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.208 | +1.009 | +0.14 pp | 1.46 pp | 2.77 pp | 2/6 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.547 | +1.885 | +0.85 pp | 1.93 pp | 6.03 pp | 3/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.773 | +3.239 | +1.94 pp | 3.62 pp | 12.31 pp | 1/5 |
| mixing_chi_0.10_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.111 | -2.076 | -0.66 pp | 1.44 pp | 2.33 pp | 1/4 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.986 | -1.062 | +0.79 pp | 0.83 pp | 2.48 pp | 3/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | +0.533 | +1.751 | -9.90 pp | 9.90 pp | 12.31 pp | 0/6 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +1.041 | +2.378 | -9.67 pp | 9.67 pp | 11.77 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +1.120 | +3.586 | -8.95 pp | 9.71 pp | 13.04 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.400 | -1.787 | -11.66 pp | 11.66 pp | 13.04 pp | 0/4 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +1.204 | -0.845 | -10.49 pp | 10.49 pp | 11.21 pp | 0/5 |
| without_skin_friction | all_digitized_points | 10000 | 6 | -1.218 | +1.198 | +2.416 | -1.58 pp | 2.97 pp | 4.31 pp | 0/6 |
| without_skin_friction | all_digitized_points | 12000 | 5 | -1.338 | +1.480 | +2.818 | -1.65 pp | 3.53 pp | 4.71 pp | 0/5 |
| without_skin_friction | all_digitized_points | 14000 | 5 | -2.466 | +1.436 | +3.902 | -1.04 pp | 5.12 pp | 10.19 pp | 0/5 |
| without_skin_friction | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.700 | -1.488 | -3.85 pp | 3.85 pp | 4.97 pp | 0/4 |
| without_skin_friction | all_digitized_points | 16000 | 5 | +2.049 | +1.437 | -0.611 | -2.71 pp | 2.71 pp | 3.43 pp | 0/5 |

## Adjacent-flow slopes for the baseline model

These secant slopes compare each neighboring pair of digitized flows. They show where the residual changes along a speed line and avoid compressing a local kink into one whole-line regression.

| rpm | Flow interval (kg/s) | Measured slope | Model slope | Slope error |
|---:|---:|---:|---:|---:|
| 10000 | 2.340–2.610 | +6.667 | +2.722 | -3.945 |
| 10000 | 2.610–3.110 | +4.000 | +1.850 | -2.150 |
| 10000 | 3.110–3.710 | -0.500 | +0.646 | +1.146 |
| 10000 | 3.710–4.280 | -8.421 | -0.712 | +7.709 |
| 10000 | 4.280–4.910 | -0.952 | -2.318 | -1.366 |
| 12000 | 2.860–3.420 | +5.179 | +2.050 | -3.129 |
| 12000 | 3.420–3.990 | +1.404 | +1.191 | -0.213 |
| 12000 | 3.990–4.540 | -0.909 | +0.381 | +1.290 |
| 12000 | 4.540–5.110 | -12.456 | -0.483 | +11.973 |
| 14000 | 3.300–3.990 | +4.928 | +1.839 | -3.088 |
| 14000 | 3.990–4.560 | +2.105 | +1.166 | -0.939 |
| 14000 | 4.560–5.110 | +2.364 | +0.588 | -1.776 |
| 14000 | 5.110–5.690 | -25.517 | -0.020 | +25.498 |
| 16000 | 3.760–4.310 | +4.182 | +1.774 | -2.408 |
| 16000 | 4.310–4.880 | +2.982 | +1.313 | -1.670 |
| 16000 | 4.880–5.440 | +0.357 | +0.872 | +0.515 |
| 16000 | 5.440–6.000 | +1.071 | +0.434 | -0.637 |

The skin-friction and chi rows are one-factor sensitivity cases; they do not establish the physical correctness of either setting.
Agreement on approximate digitized points is a model screen, not independent validation or a calibrated loss law.

Complete summary: `speedline_slope_summary.csv`; point residuals: `pointwise_efficiency_residuals.csv`; adjacent-pair slopes: `adjacent_flow_slopes.csv`.
