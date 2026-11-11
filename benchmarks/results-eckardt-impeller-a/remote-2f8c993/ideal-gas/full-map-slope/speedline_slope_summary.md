# Eckardt A efficiency slope screen across speed lines

Digitized efficiency points: 21. The plot symbols are approximate (about ±0.01 absolute efficiency).
Ablation rows are checked against the current digitization before analysis.
Slopes are unweighted least-squares fits of efficiency versus mass flow.
The full 14,000 rpm line is the primary result; a second row omits its 5.69 kg/s symbol only to show that point's influence.

| Scenario | Scope | rpm | n | Measured slope | Model slope | Slope error | Efficiency bias | MAE | Max abs | Within ±1 pp |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | all_digitized_points | 10000 | 6 | -1.218 | +0.138 | +1.355 | -4.35 pp | 4.35 pp | 6.55 pp | 0/6 |
| all_losses_enthalpy_equivalent | all_digitized_points | 12000 | 5 | -1.338 | +0.781 | +2.119 | -3.86 pp | 4.51 pp | 5.91 pp | 0/5 |
| all_losses_enthalpy_equivalent | all_digitized_points | 14000 | 5 | -2.466 | +0.940 | +3.406 | -2.94 pp | 6.00 pp | 7.66 pp | 0/5 |
| all_losses_enthalpy_equivalent | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.253 | -1.935 | -5.59 pp | 5.59 pp | 7.12 pp | 0/4 |
| all_losses_enthalpy_equivalent | all_digitized_points | 16000 | 5 | +2.049 | +1.094 | -0.955 | -4.28 pp | 4.28 pp | 5.04 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.499 | +0.719 | +3.48 pp | 3.48 pp | 5.92 pp | 0/6 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.347 | +1.685 | +4.37 pp | 4.37 pp | 9.30 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.628 | +3.095 | +5.60 pp | 5.60 pp | 15.78 pp | 0/5 |
| mixing_chi_0.05_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +0.987 | -2.201 | +3.05 pp | 3.05 pp | 5.39 pp | 0/4 |
| mixing_chi_0.05_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.890 | -1.159 | +4.60 pp | 4.60 pp | 6.38 pp | 0/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | -0.215 | +1.003 | +0.14 pp | 1.45 pp | 2.77 pp | 2/6 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +0.543 | +1.881 | +0.85 pp | 1.92 pp | 6.03 pp | 3/5 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +0.770 | +3.236 | +1.94 pp | 3.62 pp | 12.31 pp | 1/5 |
| mixing_chi_0.10_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.109 | -2.078 | -0.65 pp | 1.44 pp | 2.33 pp | 1/4 |
| mixing_chi_0.10_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +0.984 | -1.064 | +0.79 pp | 0.84 pp | 2.48 pp | 3/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 10000 | 6 | -1.218 | +0.528 | +1.746 | -9.90 pp | 9.90 pp | 12.31 pp | 0/6 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 12000 | 5 | -1.338 | +1.038 | +2.376 | -9.67 pp | 9.67 pp | 11.76 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 14000 | 5 | -2.466 | +1.118 | +3.584 | -8.94 pp | 9.70 pp | 13.04 pp | 0/5 |
| mixing_chi_0.20_bstar_1.00 | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.399 | -1.789 | -11.65 pp | 11.65 pp | 13.04 pp | 0/4 |
| mixing_chi_0.20_bstar_1.00 | all_digitized_points | 16000 | 5 | +2.049 | +1.203 | -0.846 | -10.49 pp | 10.49 pp | 11.20 pp | 0/5 |
| without_skin_friction | all_digitized_points | 10000 | 6 | -1.218 | +1.194 | +2.412 | -1.58 pp | 2.97 pp | 4.31 pp | 0/6 |
| without_skin_friction | all_digitized_points | 12000 | 5 | -1.338 | +1.478 | +2.816 | -1.65 pp | 3.53 pp | 4.71 pp | 0/5 |
| without_skin_friction | all_digitized_points | 14000 | 5 | -2.466 | +1.434 | +3.900 | -1.04 pp | 5.12 pp | 10.19 pp | 0/5 |
| without_skin_friction | 14k_sensitivity_without_5.69kg_s | 14000 | 4 | +3.188 | +1.699 | -1.489 | -3.85 pp | 3.85 pp | 4.97 pp | 0/4 |
| without_skin_friction | all_digitized_points | 16000 | 5 | +2.049 | +1.436 | -0.612 | -2.71 pp | 2.71 pp | 3.43 pp | 0/5 |

## Adjacent-flow slopes for the baseline model

These secant slopes compare each neighboring pair of digitized flows. They show where the residual changes along a speed line and avoid compressing a local kink into one whole-line regression.

| rpm | Flow interval (kg/s) | Measured slope | Model slope | Slope error |
|---:|---:|---:|---:|---:|
| 10000 | 2.340–2.610 | +6.667 | +2.720 | -3.947 |
| 10000 | 2.610–3.110 | +4.000 | +1.847 | -2.153 |
| 10000 | 3.110–3.710 | -0.500 | +0.642 | +1.142 |
| 10000 | 3.710–4.280 | -8.421 | -0.719 | +7.702 |
| 10000 | 4.280–4.910 | -0.952 | -2.329 | -1.377 |
| 12000 | 2.860–3.420 | +5.179 | +2.048 | -3.130 |
| 12000 | 3.420–3.990 | +1.404 | +1.188 | -0.215 |
| 12000 | 3.990–4.540 | -0.909 | +0.377 | +1.286 |
| 12000 | 4.540–5.110 | -12.456 | -0.489 | +11.967 |
| 14000 | 3.300–3.990 | +4.928 | +1.838 | -3.089 |
| 14000 | 3.990–4.560 | +2.105 | +1.164 | -0.941 |
| 14000 | 4.560–5.110 | +2.364 | +0.585 | -1.779 |
| 14000 | 5.110–5.690 | -25.517 | -0.023 | +25.494 |
| 16000 | 3.760–4.310 | +4.182 | +1.773 | -2.409 |
| 16000 | 4.310–4.880 | +2.982 | +1.312 | -1.671 |
| 16000 | 4.880–5.440 | +0.357 | +0.870 | +0.513 |
| 16000 | 5.440–6.000 | +1.071 | +0.432 | -0.640 |

The skin-friction and chi rows are one-factor sensitivity cases; they do not establish the physical correctness of either setting.
Agreement on approximate digitized points is a model screen, not independent validation or a calibrated loss law.

Complete summary: `speedline_slope_summary.csv`; point residuals: `pointwise_efficiency_residuals.csv`; adjacent-pair slopes: `adjacent_flow_slopes.csv`.
