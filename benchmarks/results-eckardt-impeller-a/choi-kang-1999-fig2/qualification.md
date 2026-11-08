# Eckardt A / Choi–Kang Figure 2 comparison

Three approximate experimental symbols digitized from Choi and Kang (1999), Figure 2. The test is at 14,000 rpm; mass flow is normalized by the reported design flow 4.54 kg/s. The source uses total pressure at r/r2=1.075 and total temperature at r/r2=1.687 to form its efficiency. Predictions use those same two model stations.

Digitization uncertainty is estimated at about ±0.015 in pressure ratio, temperature-rise coefficient, and efficiency. These three points are a component/performance screen, not broad map qualification.

| Flow / design flow | Flow (kg/s) | Measured PR | Predicted PR | PR error | Measured temperature coefficient | Predicted coefficient | Measured efficiency | Predicted efficiency | Efficiency error |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.75 | 3.405 | 2.010 | 1.9525 | -2.86% | 0.805 | 0.8715 | 0.920 | 0.8140 | -10.60 pp |
| 1.00 | 4.540 | 1.960 | 1.8913 | -3.51% | 0.765 | 0.8192 | 0.930 | 0.8209 | -10.91 pp |
| 1.50 | 6.810 | 1.680 | 1.6955 | +0.92% | 0.655 | 0.6935 | 0.805 | 0.7906 | -1.44 pp |

Across three points: pressure-ratio MAPE 2.43%; mean absolute temperature-coefficient error 0.0531; efficiency MAE 7.65 percentage points.

The source averages experimental quantities from different radial planes; the model comparison follows that convention. The data are approximate plot digitizations, and the report keeps their uncertainty visible rather than treating them as exact points.

Source: Y. S. Choi and S. H. Kang, ‘Application of Through-Flow Calculation to Design and Performance Prediction of Centrifugal Compressor,’ International Journal of Rotating Machinery, 5(1), 17–33 (1999), Fig. 2, DOI [10.1155/S1023621X99000020](https://doi.org/10.1155/S1023621X99000020).

Pointwise results: `qualification.csv`.

## Loss-term sensitivity at matched measurement planes

Each scenario recomputes pressure at r/r2=1.075 and total temperature at r/r2=1.6875, then forms mixed-plane efficiency using those same stations. Rows are sensitivity experiments, not evidence that removing a loss term is physically correct.

| Scenario | Points | Efficiency bias (pp) | Efficiency MAE (pp) | Max abs (pp) | Mean ΔT (K) | Mean near PR |
|---|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | 3 | -7.65 | 7.65 | 10.91 | 68.02 | 1.8464 |
| shroud_inlet_speed_proxy | 3 | -9.17 | 9.17 | 12.35 | 67.79 | 1.8239 |
| legacy_mean_temperature_entropy | 3 | -8.68 | 8.68 | 11.71 | 67.85 | 1.8314 |
| rectangular_endpoint_hydraulic_diameter | 3 | -6.46 | 6.74 | 9.95 | 68.22 | 1.8639 |
| without_skin_friction | 3 | -2.65 | 6.77 | 7.34 | 68.84 | 1.9215 |
| without_blade_loading | 3 | -6.74 | 6.74 | 9.91 | 68.13 | 1.8608 |
| without_incidence | 3 | -6.94 | 7.18 | 10.91 | 68.16 | 1.8565 |
| without_clearance | 3 | -7.33 | 7.33 | 10.60 | 68.07 | 1.8513 |
| without_mixing | 3 | +2.78 | 2.81 | 7.58 | 60.82 | 1.8573 |
| without_disk_friction | 3 | -7.36 | 7.36 | 10.62 | 67.80 | 1.8468 |
| without_mixing_and_disk_friction | 3 | +3.14 | 3.14 | 7.83 | 60.60 | 1.8576 |
| without_impeller_losses | 3 | -0.71 | 6.38 | 8.52 | 69.12 | 1.9521 |
| without_diffuser_friction | 3 | -7.32 | 7.32 | 10.59 | 68.02 | 1.8508 |
| mixing_chi_0.05_bstar_1.00 | 3 | +0.70 | 3.16 | 5.80 | 62.13 | 1.8553 |
| mixing_chi_0.10_bstar_1.00 | 3 | -2.87 | 4.68 | 5.94 | 64.52 | 1.8516 |
| mixing_chi_0.20_bstar_1.00 | 3 | -13.54 | 13.54 | 17.03 | 72.91 | 1.8394 |
| mixing_chi_0.25_bstar_1.00 | 3 | -20.41 | 20.41 | 24.31 | 79.61 | 1.8303 |
| mixing_chi_0.30_bstar_1.00 | 3 | -28.05 | 28.05 | 32.40 | 88.73 | 1.8186 |
| mixing_chi_0.15_bstar_1.02 | 3 | -8.68 | 8.68 | 11.99 | 68.83 | 1.8453 |
| mixing_chi_0.15_bstar_1.10 | 3 | -13.07 | 13.07 | 16.55 | 72.50 | 1.8400 |
| correlation_exit_angle_60deg_screen | 3 | -10.46 | 10.46 | 13.12 | 67.54 | 1.8061 |

Pointwise ablation: `station_matched_loss_ablation.csv`; summary: `station_matched_loss_ablation_summary.csv`.
