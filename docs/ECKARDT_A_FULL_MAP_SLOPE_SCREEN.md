# Eckardt A efficiency slope screen

The 14,000 rpm finite-difference diagnostic suggested a slope mismatch between
the model and the digitized efficiency symbols. This screen checks whether the
mismatch also appears on the other digitized speed lines and whether the
skin-friction and mixing-chi sensitivities have similar effects there.

The script consumes the benchmark's `loss_term_ablation.csv`, verifies that its
measured efficiency values match the current digitized map, and reports a
per-speed-line measured/model efficiency slope, slope residual, bias, MAE, and
maximum absolute error. It includes the baseline, no-skin-friction, and
available mixing-chi scenarios. It preserves all digitized points as the
primary result and reports a separate 14,000 rpm sensitivity row omitting the
5.69 kg/s point.

The measured symbols are approximate plot digitizations (about ±0.01 absolute
efficiency). The slope calculations are unweighted. Scenario changes are
sensitivity experiments, not evidence that a loss term should be removed or
that a parameter setting is physically correct.

The source labels the plotted isentropic efficiency as total-to-total
efficiency (ηtt). Its design-validation table reports Impeller A experimental
efficiency 0.876 at 4.54 kg/s; the digitized 14,000 rpm value 0.877 at 4.56
kg/s is consistent with that tabulated point. This supports the efficiency
definition and the near-design digitization. The 5.69 kg/s symbol is not
tabulated, so its value remains a plot read and should be checked against the
figure before treating its local residual as a physical feature. The paper's
table and figure identify the A point and efficiency metric, but do not
independently validate this implementation. (Meroni et al., *Applied Energy*
232 (2018), 139–156, DOI: 10.1016/j.apenergy.2018.09.210.)

The analysis also writes every point residual and every adjacent-flow secant
slope. Those outputs localize rapid residual changes that a single speed-line
regression can hide. Adjacent slopes remain sensitive to digitization noise,
especially where flow points are close together.

Generate fresh ablation data using the corrected efficiency file, then run:

```bash
python3 scripts/analyze_eckardt_full_map_slope.py \
  --ablation-csv benchmarks/results-eckardt-impeller-a/corrected-digitization/loss_term_ablation.csv \
  --efficiency-map data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  --output-dir benchmarks/results-eckardt-impeller-a/full-map-slope
```

The report and machine-readable outputs are written to
`speedline_slope_summary.md`, `speedline_slope_summary.csv`,
`pointwise_efficiency_residuals.csv`, and `adjacent_flow_slopes.csv` in the
output directory. A stale ablation file with the old 16,000 rpm, 6.00 kg/s
value is rejected instead of being silently compared with the corrected map.
