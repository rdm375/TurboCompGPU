# Eckardt A impeller-loss screening

## Purpose

This experiment adds a first physics-based impeller loss estimate to the
Eckardt A mean-line benchmark. It checks pressure ratio and total-to-total
efficiency together, while retaining the smooth-wall diffuser estimate as the
reference case. It is a screening model, not compressor-map qualification.

The pressure-ratio map and efficiency symbols were digitized independently.
They are not pointwise paired measurements, and the efficiency values have
approximate digitization uncertainty of 0.01 absolute. In particular, the
high-flow symbols at 14,000 rpm / 5.69 kg/s and 16,000 rpm / 6.00 kg/s were
visually checked against Figure 6(d) of Meroni et al.; they lie near ηtt=0.742
and ηtt=0.655 respectively. The abrupt high-flow efficiency fall remains a
large model discrepancy.

## Model added

The benchmark estimates two impeller losses from the operating point:

- **Skin friction:** Jansen/Aungier form, using the Jansen blade-path estimate
  and a Darcy friction factor. The channel hydraulic diameter is estimated
  from rectangular inlet and exit passages built from the published annulus,
  blade pitch, width, and blockage. The friction factor is a smooth-wall
  Petukhov estimate at the passage Reynolds number.
- **Blade loading:** Coppage correlation
  `Δh_bl = 0.05 D_f² U₂²`, with `C_df = 0.75` for the un-splittered Eckardt A
  impeller. The diffusion factor uses the published geometry and the operating
  Euler work.

The sum of these losses is represented as entropy generation
`Δs = (Δh_sf + Δh_bl) / T_mean` and passed through the existing operating-point
solver. The loss estimate is recalculated from the corrected state until the
pressure change converges. Diffuser friction remains independent: the smooth
wall Reynolds-number estimate is the reference; fixed Darcy factors from 0 to
0.04 are sensitivity cases only.

The equations follow the loss-correlation set summarized by Meroni et al.,
Table 4, including the Jansen path-length approximation, Coppage blade-loading
loss, and the listed Eckardt A dimensions. The implementation's rectangular
hydraulic-diameter approximation is an explicit geometry simplification.
Incidence, tip-clearance leakage, jet-wake mixing, disk friction, recirculation,
and shock losses are not included. No factor was fitted to these digitized
data.

## Results

Using the smooth-wall diffuser reference:

| Metric | No impeller-loss terms | With skin-friction and blade-loading terms |
|---|---:|---:|
| Pressure-ratio MAPE (17 points) | 2.39% | 2.23% |
| Pressure-ratio maximum absolute error | 8.78% | 3.47% |
| Efficiency mean signed error (21 points) | +12.23 pp | +5.18 pp |
| Efficiency MAE (21 points) | 12.23 pp | 5.18 pp |
| Efficiency maximum absolute error | 31.10 pp | 23.98 pp |

The model improves both aggregate comparisons, particularly the efficiency
MAE, without tuning the diffuser friction factor. It still overpredicts
efficiency by about 24 percentage points at 16,000 rpm and 6.00 kg/s. This
remaining high-flow error is consistent with important omitted off-design
physics, but the current data and model do not identify which omitted
mechanism dominates. The fixed-friction sensitivity sweep also shows that
selecting friction to reduce efficiency error can worsen pressure-ratio error;
those cases must not be treated as calibration.

The pointwise output CSVs now include predicted ηtt, impeller skin-friction and
blade-loading losses, the skin-friction factor, and the diffusion factor.

## Reproduce

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  benchmarks/results-eckardt-impeller-a/impeller_loss_pressure_ratio.csv \
  data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  benchmarks/results-eckardt-impeller-a/impeller_loss_efficiency.csv
```

## Source

Meroni, A., et al. “Design of centrifugal compressors for heat pump systems.”
*Applied Energy* 232 (2018), 657–667. The paper's Table 4 summarizes the
impeller loss correlations and Table 5 gives the Eckardt A geometry and
clearances. [DTU repository postprint (PDF)](https://backend.orbit.dtu.dk/ws/files/154611324/Postprint_Design_of_centrifugal_compressors_for_heat_pump_systems.pdf).
