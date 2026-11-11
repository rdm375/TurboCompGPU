# Published Eckardt compressor loss-set diagnostics

## Purpose

The Eckardt A friction-sweep benchmark now evaluates two correlation bundles
reported by Kovář et al. (2021): Set 2042, selected as their best set for
Eckardt type A, and Set 2137, their cross-geometry compromise. They are
additional diagnostic scenarios in `loss_term_ablation.csv`; the benchmark's
current screen remains the default. These runs do not reproduce the paper's
full prediction algorithm or qualify either set for production use.

Kovář et al. compared 4,608 combinations of individual correlations. Their
paper identifies type A at 14,000 and 16,000 rpm as difficult off-design
conditions and describes a reconstructed three-dimensional geometry using
blade-thickness distributions and heel radii. The current benchmark is a
one-dimensional reconstruction with materially less geometry detail.

## Correlations in the two sets

| Loss | Set 2042 (Eckardt A) | Set 2137 (cross-geometry) |
|---|---|---|
| Incidence | Conrad | Conrad |
| Choke | Aungier | Aungier |
| Entrance diffusion | Aungier | Aungier |
| Blade loading | Aungier | Aungier |
| Skin friction | Jansen | Aungier velocity average in Jansen form |
| Clearance | Rodgers | Jansen |
| Mixing | Aungier | Aungier |
| Recirculation | Coppage | Coppage |
| Disk friction | Shepherd | Shepherd |
| Vaneless diffuser | Coppage | Stanitz |

The correlation definitions follow Kovář et al., Sections 3.1–3.3 and Table 3.
In particular, Aungier's blade loading is `ΔW²/48`; Coppage recirculation is
`0.02 sqrt(tan(alpha₂)) D_f² U₂²`; Shepherd disk friction uses
`0.01356 rho₂ U₂³ D₂²/(m_dot Re₂^0.2)`; and Rodgers clearance is
`0.1 (epsilon/b₂) U₂²`.

## Implementation and limits

- Conrad incidence uses `sigma=0.6` and the reported type-A meanline inlet angle (56°), inlet radius, and blade
  thickness to estimate the zero-incidence metal angle. The incidence
  triangle remains the benchmark's meanline triangle.
- The Aungier throat-shock loss uses the paper's polynomial. The exact
  `A_th*/(C_r A_th)` ratio is not available from the current coupled throat
  solver, so the benchmark substitutes inverse throat-capacity utilization.
  This is explicitly a proxy; it is not the paper's corrected sonic-area
  calculation.
- Aungier entrance diffusion uses the coupled operating-throat relative speed
  for `W_th`. The paper's method solves a separate three-equation throat model;
  the existing throat state is an approximation to that input.
- Jansen skin friction uses the existing Schlichting Fanning factor and
  meanline hydraulic diameter/path length. Set 2137 changes the averaged
  velocity to Aungier's RMS expression.
- The source type-A clearance input is 0.5 mm. Rodgers and Jansen use that
  clearance; the baseline screen retains its earlier separate clearance
  assumptions.
- Aungier mixing follows Equations (43)–(45): `Wmax` uses the mean-streamline
  relative speed `W1xi`, and `Wout` uses the inlet relative tangential component
  `Wu1xi`. This keeps those inputs distinct from the shroud speed `W1s` used
  by other loss correlations.
- The two diffuser correlations return loss rather than the state trajectory
  used by this repository's variable-width diffuser solver. The benchmark
  matches each correlation's outlet entropy loss with an equivalent Darcy
  factor and then re-integrates its mass, energy, and entropy equations. The
  reported `darcy_factor` for these scenarios is therefore an equivalent
  numerical parameter, not a measured wall-friction factor. Coppage's target
  uses Equation (58); Stanitz's target is evaluated on the smooth-wall reference
  diffuser state before being mapped to an equivalent factor.
- The Stanitz correlation is evaluated with `cp=1004.5 J/(kg K)` and the
  corresponding ideal-air `gamma`; the resulting entropy target is applied to
  whichever property model (`ideal-gas` or `pr-air`) the run selected.
- The current diffuser outlet taper, inlet/outlet radii, blockage, and
  velocity solution remain those in the benchmark. They are not silently
  replaced with the paper's full geometry reconstruction.

The new scenarios report the Aungier choke, entrance-diffusion, and Coppage
recirculation terms in the appended `impeller_choke_loss_j_kg`,
`impeller_entrance_diffusion_loss_j_kg`, and
`impeller_recirculation_loss_j_kg` columns. Inspect those values alongside
pressure ratio, efficiency, and the component losses. Do not infer validation
from an aggregate score or from a closer match to the same digitized points.

## Reproduce

```sh
cmake --build build-model-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-model-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  benchmarks/results-eckardt-impeller-a/published-loss-sets/pressure.csv \
  data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  benchmarks/results-eckardt-impeller-a/published-loss-sets/efficiency.csv \
  pr-air
```

The benchmark writes the two set rows to the sibling
`loss_term_ablation.csv` and preserves the selected property model in
`thermodynamic_model.txt`.

## Source

P. Kovář, A. Tater, P. Mačák, and T. Vampola, “Searching for the Most Suitable
Loss Model Set for Subsonic Centrifugal Compressors Using an Improved Method
for Off-Design Performance Prediction,” *Energies* 14(24), 8545 (2021),
[doi:10.3390/en14248545](https://doi.org/10.3390/en14248545). Set membership is
from Table 3; the Galvas pressure update and the internal/external loss split
are described in Section 2; the selected correlations are in Sections 3.1–3.3.

The publisher PDF is archived at
[`docs/references/kovar-tater-macak-vampola-2021-eckardt-loss-sets.pdf`](references/kovar-tater-macak-vampola-2021-eckardt-loss-sets.pdf).
It is reproduced under the article's Creative Commons Attribution 4.0 license.
The paper's Set 2042 is a set of correlations evaluated inside the authors'
Galvas-based algorithm and reconstructed geometry. This benchmark implements
the selected correlation bundle with real-gas state calculations and maps the
diffuser target to its existing one-dimensional solver; it is not a reproduction
of the paper's complete algorithm or geometry reconstruction.
