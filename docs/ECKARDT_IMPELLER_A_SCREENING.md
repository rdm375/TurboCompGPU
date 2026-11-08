# Eckardt Impeller A map screening comparison

## Purpose and status

This is the first external-data comparison for the complete-stage solver. It
uses 17 Eckardt Impeller A map points at 10,000, 12,000, 14,000, and 16,000
rpm. The tabulated flow/pressure-ratio values are reproduced from Table 5 of
Bayomi et al., which labels them experimental data. The original experimental
work is Eckardt's compressor study, also identified by the public validation
implementation as Meroni et al. (2018), DOI
[10.1016/j.apenergy.2018.09.210](https://doi.org/10.1016/j.apenergy.2018.09.210).

This is a **screening comparison**, not a qualification claim. Bayomi et al.'s
map table labels the output as pressure ratio without subscripts. Meroni et
al.'s primary validation paper labels the corresponding map quantity
`PRtt`, total-to-total pressure ratio, so that is the comparison used here.

> **Superseded results:** the metrics below predate correction of the inlet
> annulus radii to the diameter fields expected by the model. Use the corrected
> geometry results and friction sensitivity in
> [ECKARDT_DIFFUSER_FRICTION_SCREENING.md](ECKARDT_DIFFUSER_FRICTION_SCREENING.md).

## Reproduction

Configure with a ThermoGPU installation, then build and run:

```sh
cmake -S . -B build-eckardt-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/ThermoGPU/install

cmake --build build-eckardt-validation \
  --target turbocompgpu_validate_eckardt_impeller_a --parallel 8

./build-eckardt-validation/turbocompgpu_validate_eckardt_impeller_a
```

The executable reads `data/validation/eckardt_impeller_a.csv` and writes
`benchmarks/results-eckardt-impeller-a/qualification.csv`. It reports failures
as rows too, so failed points are not silently removed from the record.
An optional third argument overrides the slip factor for controlled
sensitivity checks, for example:

```sh
./build-eckardt-validation/turbocompgpu_validate_eckardt_impeller_a \
  data/validation/eckardt_impeller_a.csv /tmp/eckardt-sigma1.csv 1.0
```

## Baseline result

Run with GCC 15.2 and the CPU ThermoGPU 2.x installation on 2026-11-07. The
first run used `r3/r2=1.69`, copied from a secondary implementation. After
correcting it to the primary paper's Impeller A geometry (`r3/r2=2.69`), the
reported total-to-total pressure ratios were unchanged at the displayed
precision, as expected for this zero-loss diffuser model. The per-point CSV
was regenerated using the primary geometry.

| Speed | Points solved | Signed bias | Mean absolute relative error | Maximum absolute error |
|---:|---:|---:|---:|---:|
| 10,000 rpm | 5/5 | +0.88% | 0.98% | 2.01% |
| 12,000 rpm | 4/4 | +0.01% | 1.46% | 2.89% |
| 14,000 rpm | 4/4 | −3.23% | 5.03% | 7.22% |
| 16,000 rpm | 4/4 | −6.96% | 7.15% | 11.13% |
| **All points** | **17/17** | **−2.14%** | **3.50%** | **11.13%** |

No parameters were fitted against these points. The residual is systematically
negative on the higher speed lines: the maximum error is −11.13% at 16,000 rpm
and 6.8 kg/s. For context, Meroni et al.'s more detailed model reports up to
6% pressure-ratio underprediction for Impeller A and notes that it does not
reproduce the speed-line trends perfectly.

## Model and assumptions

- Air is modeled with ThermoGPU's default ideal gas (`R=287.05 J/(kg K)`,
  `cp=1004.5 J/(kg K)`). The source reports inlet total conditions
  `T01=288.15 K`, `p01=1.01 bar`; the harness converts these to a static inlet
  state and velocity on the subsonic ideal-gas branch at each mass flow. Inlet
  prewhirl is zero.
- Meroni et al.'s primary geometry table supplies 20 blades, `r2=0.2 m`,
  `r1,shroud=0.14 m`, `r1,hub=0.06 m`, and `b2=0.026 m`.
- Exit blade thickness is `1.08 mm`, per Meroni et al.; the slip-factor
  correlation is applied using 20 blades and a 30-degree exit angle.
- The exit-angle conversion to TurboCompGPU's convention is 60 degrees. The
  diffuser uses `r3/r2=2.69`, `b3/b2=0.51`, and zero blockage, per the primary
  geometry table.
- The inlet annulus geometry is represented directly by its hub and shroud
  radii. This simplified model does not represent the inducer blade throat.
- Entropy generation and diffuser blockage are zero. The current complete
  stage model therefore has no calibrated impeller, diffuser, leakage, disk,
  or volute losses.
- The comparison metric uses only converged predictions. All 17 points
  converged in this run; a future run with failures must report coverage beside
  error metrics.

The CSV also records the impeller total-pressure ratio and the diffuser ratio
`p03/p02`. That diffuser ratio is exactly 1.0 for all 17 points, confirming
that this lossless diffuser contributes no total-pressure change; the map
discrepancy is entirely present at the impeller in this model.

The per-point predictions are preserved in
[`benchmarks/results-eckardt-impeller-a/qualification.csv`](../benchmarks/results-eckardt-impeller-a/qualification.csv).

## Controlled slip sensitivity

Setting the slip factor to its upper bound, `sigma=1.0`, is a sensitivity
experiment, not a physical fit. All 17 points still converge. The aggregate
MAPE rises from 3.50% to 5.16% and signed bias changes from −2.14% to +4.12%.
High-flow errors improve (for example, 16,000 rpm and 6.8 kg/s changes from
−11.13% to −4.47%), while low-flow errors become positive and reach +12.41%.
One constant slip factor therefore cannot correct the full map.

## Exit-whirl equation comparison

Meroni et al.'s published exit-whirl relation is
`Ctheta2 = sigma*U2 - Cm2*tan(abs(beta2b))`. With the angle conversion used
here, that is `sigma*U2 - Cm2/tan(beta_cpp)`. The current TurboCompGPU
implementation instead evaluates `sigma*(U2 - Cm2/tan(beta_cpp))`. For
`sigma < 1`, the current formula predicts greater exit whirl by
`(1-sigma)*Cm2/tan(beta_cpp)`. This difference has the wrong sign to explain
the observed pressure-ratio underprediction: matching the published relation
would reduce the current model's predicted work further. The work/slip closure
must be resolved before calibrating loss terms.

These assumptions make the result useful for checking solver behavior and
quantifying the current gap. Changing diffuser radius did not change the total
pressure ratio because the implemented diffuser is lossless; it does change its
static outlet state. The increasing high-speed, high-flow underprediction
cannot be explained by adding losses alone. The exit-whirl closure and missing
inlet-throat/incidence physics need review before any calibration. This
comparison does **not** qualify the model.

## Next qualification steps

1. Resolve the exit-whirl formula convention against a standards-backed slip
   correlation and make the selected closure explicit in code and tests.
2. Obtain the raw source data if accessible and compare its point values to
   Bayomi et al.'s transcription and Meroni et al.'s plotted map.
3. Add explicit data provenance, unit checks, and an immutable source checksum.
4. Add inlet-throat geometry and incidence/inducer loss physics from the
   primary geometry inputs.
5. Compare more observables if available: temperature rise, efficiency, and
   diffuser/exit states, not pressure ratio alone.
6. Repeat against Impeller B or O without retuning the Impeller A geometry or
   loss assumptions; keep calibration and held-out datasets separate.
7. Add an acceptance policy before treating any map as qualified, including
   coverage requirements and error tolerances by speed line and operating
   region.
