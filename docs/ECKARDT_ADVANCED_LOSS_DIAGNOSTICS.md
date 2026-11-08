# Eckardt A angle, throat-speed, and mixing diagnostics

## Scope

This screening pass compares the 30° and 60° outlet-angle conventions in the
Galvas hydraulic-diameter and Jansen blade-path correlations, replaces the
shroud-inlet proxy in Jansen's mean-relative-speed estimate with an
operating-point throat solution, and varies the Johnston-Dean mixing inputs.
The compressor velocity triangle continues to use TurboCompGPU's 60° angle.

The pressure map has 17 points. The efficiency data has 21 approximate,
independently digitized points, so these are separate comparisons. No
parameter was fitted to the data.

## Findings

### Outlet-angle convention affects skin friction

The source geometry gives an outlet metal angle of -30° in its convention,
while TurboCompGPU stores the complementary 60° velocity-triangle angle. The
loss correlations that use `cos(beta2)` now use the 30° magnitude; the 60°
correlation case remains as a comparison.

Changing the correlation angle from 60° to 30° raises the mean Galvas passage
diameter from 29.69 mm to 33.05 mm and lowers the mean Jansen path from
307.08 mm to 230.93 mm. With the operating throat speed, the 30° case has
lower pressure-ratio MAPE (2.74% versus 4.44%) and lower efficiency MAE
(7.40 versus 9.45 percentage points).

### The operating throat speed reduces the friction estimate

The previous estimate used the shroud inlet relative speed `W1s` as a proxy
for Jansen's throat speed `W_th,1`. The new calculation solves a subcritical
blade-aligned throat state at each operating point. It conserves inlet
entropy and rotor rothalpy and solves `mdot = rho A_eff W_n`; the sonic state
provides the upper bracket. It uses the existing meanline throat geometry,
which does not locate the true three-dimensional minimum passage.

The throat area initially used the source's meridional angle directly in a
helper that expects angle from the tangential direction. The angles are now
converted explicitly: the mean 51.5-degree meridional angle becomes 38.5
degrees for throat-area projection. The area estimate falls from 0.0390 to
0.0303 m², and the 21-point mean throat speed rises from 82.65 to 109.03 m/s
(range 60.80–156.52 m/s). The selected mean relative speed `Wbar` stays nearly
unchanged at 140.44 m/s because the mean-inlet branch remains the larger
candidate. Mean skin-friction loss is 2.74 kJ/kg. On the same screening data,
pressure-ratio MAPE is 2.58% and efficiency MAE is 7.07 pp, with maximum
absolute efficiency error 13.71 pp.

The independent audit checks throat continuity, rothalpy, entropy, and
subsonic conditions as well as `Wbar` and skin-loss arithmetic. Across 441
scenario rows, maximum relative residuals are `4.91e-12` for throat
continuity, `1.22e-11` for rothalpy, and `7.22e-12` for `Wbar`; maximum
absolute entropy residual is `9.74e-11 J/(kg K)`. The maximum throat Mach
number is 0.457. These checks confirm the implementation's equations under
its meanline, adiabatic, isentropic assumptions; they do not validate the
meanline throat geometry against a measured blade passage.

### The mixing term is highly sensitive to its inputs

The Johnston-Dean `b*` is diffuser-inlet depth divided by impeller-exit depth.
For this geometry both physical depths are about 26 mm, so `b*` is
approximately 1. The `b3/b2=0.51` value describes the downstream diffuser
width taper; it is not Johnston-Dean `b*`. The cited sudden-expansion model
is stated for `b* > 1`, making `b*=1` a boundary case. A contraction needs a
wall-shape-dependent coefficient that this model does not include.

At `b*=1`, changing `chi` from 0.15 to 0.10 and 0.05 reduces efficiency MAE
from 7.07 pp to 3.64 and 3.30 pp. Raising it to 0.20 and 0.30 increases MAE
to 11.74 and 25.04 pp. These are sensitivities on the same digitized map,
not independent evidence for changing the model coefficient.

| Screen | Efficiency MAE | Pressure-ratio MAPE | Interpretation |
|---|---:|---:|---|
| `chi=0.15`, `b*=1`, 30° correlation angle, corrected throat angle | 7.07 pp | 2.58% | Current screening case |
| Same case with shroud-speed proxy | 8.27 pp | 3.73% | Diagnostic proxy |
| 60° correlation angle, corrected throat angle | 9.12 pp | 4.44% | Angle-convention comparison |
| `chi=0.05`, `b*=1` | 3.30 pp | 2.38% | Sensitivity only; not independently qualified |
| `chi=0.10`, `b*=1` | 3.64 pp | 2.53% | Sensitivity only |
| `chi=0.20`, `b*=1` | 11.74 pp | 3.02% | Sensitivity only |
| `chi=0.30`, `b*=1` | 25.04 pp | 3.87% | Sensitivity only |
| `chi=0.15`, `b*=1.02` | 7.88 pp | 2.79% | Just inside stated domain |
| `chi=0.15`, `b*=1.10` | 11.34 pp | 3.00% | Inside stated domain |
| No mixing | 4.80 pp | 2.30% | Ablation only; does not qualify the physics |

The current default mixing term lowers average predicted efficiency too much
on this digitized map. Lower `chi` improves average error but retains a large
worst-point error; it should not be selected by fitting these same points.
Independent wake measurements or a separately established
operating-point-dependent secondary-flow estimate are needed to qualify it.

## Reproducibility

The sweep writes `loss_term_ablation.csv` alongside the efficiency output. Its
rows include the operating throat state and residuals, both inlet-speed
estimates, `Wbar`, friction geometry, `chi`, `b*`, and component losses.

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/Sync/ThermoGPU/install-release-cpu"

cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep \
  data/validation/eckardt_impeller_a.csv \
  benchmarks/results-eckardt-impeller-a/operating_throat_pressure.csv \
  data/validation/eckardt_impeller_a_efficiency_digitized.csv \
  benchmarks/results-eckardt-impeller-a/operating_throat_efficiency.csv

python3 scripts/audit_eckardt_skin_friction.py \
  benchmarks/results-eckardt-impeller-a/loss_term_ablation.csv
python3 scripts/analyze_eckardt_loss_ablation.py
```

## Sources

- Meroni et al., “Design of centrifugal compressors for heat pump systems,”
  *Applied Energy* 232 (2018), 139–156. [DTU postprint](https://backend.orbit.dtu.dk/ws/files/154611324/Postprint_Design_of_centrifugal_compressors_for_heat_pump_systems.pdf).
- Johnston and Dean, “Losses in Vaneless Diffusers of Centrifugal Compressors
  and Pumps: Analysis, Experiment, and Design,” *Journal of Engineering for
  Gas Turbines and Power* (1966), [DOI: 10.1115/1.3678477](https://doi.org/10.1115/1.3678477).
- A published derivation defines `b*` as diffuser-inlet depth divided by
  impeller-exit depth and limits the sudden-expansion treatment to `b* > 1`:
  [SciELO article](https://www.scielo.org.za/scielo.php?pid=S2309-89882005000100001&script=sci_arttext).
