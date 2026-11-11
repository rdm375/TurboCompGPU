# Eckardt A / Choi–Kang Figure 2 comparison

Three approximate experimental symbols digitized from Choi and Kang (1999), Figure 2. The test is at 14,000 rpm; mass flow is normalized by the reported design flow 4.54 kg/s. The source uses total pressure at r/r2=1.075 and total temperature at r/r2=1.687 to form its efficiency. Predictions use those same two model stations.

Digitization uncertainty is estimated at about ±0.015 in pressure ratio, temperature-rise coefficient, and efficiency. These three points are a component/performance screen, not broad map qualification.

| Flow / design flow | Flow (kg/s) | Measured PR | Predicted PR | PR error | Measured temperature coefficient | Predicted coefficient | Measured efficiency | Predicted efficiency | Efficiency error |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.75 | 3.405 | 2.010 | 1.9806 | -1.46% | 0.805 | 0.8735 | 0.920 | 0.8312 | -8.88 pp |
| 1.00 | 4.540 | 1.960 | 1.9247 | -1.80% | 0.765 | 0.8224 | 0.930 | 0.8423 | -8.77 pp |
| 1.50 | 6.810 | 1.680 | 1.7507 | +4.21% | 0.655 | 0.7035 | 0.805 | 0.8306 | +2.56 pp |

Across three points: pressure-ratio MAPE 2.49%; mean absolute temperature-coefficient error 0.0581; efficiency MAE 6.74 percentage points.

The source averages experimental quantities from different radial planes; the model comparison follows that convention. The data are approximate plot digitizations, and the report keeps their uncertainty visible rather than treating them as exact points.

Source: Y. S. Choi and S. H. Kang, ‘Application of Through-Flow Calculation to Design and Performance Prediction of Centrifugal Compressor,’ International Journal of Rotating Machinery, 5(1), 17–33 (1999), Fig. 2, DOI [10.1155/S1023621X99000020](https://doi.org/10.1155/S1023621X99000020).

Pointwise results: `qualification.csv`.

## Loss-term sensitivity at matched measurement planes

Each scenario recomputes pressure at r/r2=1.075 and total temperature at r/r2=1.6875, then forms mixed-plane efficiency using those same stations. Rows are sensitivity experiments, not evidence that removing a loss term is physically correct.

| Scenario | Points | Efficiency bias (pp) | Efficiency MAE (pp) | Max abs (pp) | Mean ΔT (K) | Mean near PR |
|---|---:|---:|---:|---:|---:|---:|
| all_losses_enthalpy_equivalent | 3 | -5.03 | 6.74 | 8.88 | 68.45 | 1.8853 |
| published_eckardt_a_set_2042 | 3 | +4.09 | 5.97 | 14.73 | 62.05 | 1.8878 |
| published_cross_geometry_set_2137 | 3 | +1.48 | 6.07 | 11.32 | 61.77 | 1.8530 |
| distributed_skin_friction_screen | 3 | -5.05 | 6.75 | 8.90 | 68.45 | 1.8851 |
| geometry_continuity_skin_friction_screen | 3 | -4.51 | 6.79 | 8.59 | 68.54 | 1.8931 |
| shroud_inlet_speed_proxy | 3 | -5.91 | 7.01 | 9.75 | 68.32 | 1.8718 |
| legacy_mean_temperature_entropy | 3 | -5.67 | 6.79 | 9.42 | 68.35 | 1.8757 |
| rectangular_endpoint_hydraulic_diameter | 3 | -4.50 | 6.76 | 8.55 | 68.54 | 1.8933 |
| without_skin_friction | 3 | -2.75 | 6.82 | 7.46 | 68.82 | 1.9200 |
| without_blade_loading | 3 | -3.98 | 5.92 | 7.63 | 68.58 | 1.9020 |
| without_incidence | 3 | -4.35 | 7.20 | 8.77 | 68.58 | 1.8952 |
| without_clearance | 3 | -4.71 | 6.67 | 8.61 | 68.50 | 1.8902 |
| without_mixing | 3 | +5.77 | 5.77 | 12.14 | 61.18 | 1.8966 |
| without_disk_friction | 3 | -4.73 | 6.59 | 8.50 | 68.23 | 1.8857 |
| without_mixing_and_disk_friction | 3 | +6.14 | 6.14 | 12.40 | 60.96 | 1.8970 |
| without_impeller_losses | 3 | -0.71 | 6.38 | 8.52 | 69.12 | 1.9521 |
| without_diffuser_friction | 3 | -4.69 | 6.57 | 8.45 | 68.45 | 1.8899 |
| audit_blade_loading_shroud_speed | 3 | -5.64 | 7.00 | 9.52 | 68.37 | 1.8759 |
| audit_blade_loading_06_shroud_speed | 3 | -5.53 | 6.96 | 9.41 | 68.38 | 1.8776 |
| audit_incidence_mean_speed | 3 | -4.81 | 6.87 | 8.77 | 68.49 | 1.8885 |
| audit_clearance_exit_meridional_speed | 3 | -5.06 | 6.74 | 8.90 | 68.45 | 1.8849 |
| audit_disk_diameter_gap_ratio | 3 | -4.98 | 6.71 | 8.81 | 68.41 | 1.8854 |
| audit_combined_alternative_forms | 3 | -5.29 | 7.06 | 9.29 | 68.38 | 1.8804 |
| mixing_chi_0.05_bstar_1.00 | 3 | +3.62 | 3.62 | 10.25 | 62.50 | 1.8945 |
| mixing_chi_0.10_bstar_1.00 | 3 | -0.09 | 4.73 | 6.97 | 64.92 | 1.8907 |
| mixing_chi_0.20_bstar_1.00 | 3 | -11.12 | 11.12 | 15.37 | 73.39 | 1.8781 |
| mixing_chi_0.25_bstar_1.00 | 3 | -18.22 | 18.22 | 22.88 | 80.16 | 1.8686 |
| mixing_chi_0.30_bstar_1.00 | 3 | -26.13 | 26.13 | 31.14 | 89.37 | 1.8564 |
| mixing_chi_0.15_bstar_1.02 | 3 | -6.10 | 7.17 | 10.02 | 69.27 | 1.8841 |
| mixing_chi_0.15_bstar_1.10 | 3 | -10.64 | 10.64 | 14.86 | 72.97 | 1.8787 |
| correlation_exit_angle_60deg_screen | 3 | -6.25 | 6.66 | 9.76 | 68.25 | 1.8672 |

Pointwise ablation: `station_matched_loss_ablation.csv`; summary: `station_matched_loss_ablation_summary.csv`.

## Mixing flow diagnostic

Not run: the current benchmark source does not emit the dense chi sweep required by the historical flow-dependence diagnostic. The measured-plane comparison and the emitted named loss scenarios above are complete.
