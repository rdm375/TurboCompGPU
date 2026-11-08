# Eckardt radial-discharge Impeller O traverse: cross-family screen

## Dataset identity and comparison

`data/validation/eckardt_radial_impeller_o_diffuser_traverse.csv` transcribes the six 14,000 rpm integrated traverse points in Table 4 of D. Japikse, *Assessment of Single- and Two-Zone Modeling of Centrifugal Compressors*, ASME 85-GT-73 (1985), DOI: [10.1115/85-GT-73](https://doi.org/10.1115/85-GT-73).

**Not Impeller A data.** Japikse identifies this detailed 14,000 rpm dataset as radial-discharge Impeller O and says comparable public traverse data for the backswept impeller was unavailable. The benchmark uses backswept Impeller A geometry, so this is only a cross-family screen; its residuals cannot qualify or diagnose Impeller A.

Table 4 reports mass-averaged total pressure and flow angle at station 4 (`r/r2 = 1.075`), and static and total pressure at station 5 (`r/r2 = 1.6875`). It gives both orifice mass flow and flow integrated from the traverse. The comparator evaluates the model at the orifice mass flow, the operating flow measurement, and preserves the traverse-integrated flow in the pointwise output for context. The model angle is measured from the meridional direction.

Japikse states that the detailed traverse profiles were supplied directly by Eckardt and covers a single 14,000 rpm speed line. This public table provides integrated averages, not the full spanwise profiles or a multi-speed traverse dataset. This is a component-level cross-check of meanline predictions, not a complete compressor qualification.

## Reproduce

```sh
cmake --build build-eckardt-validation \
  --target turbocompgpu_eckardt_diffuser_friction_sweep --parallel 8

python3 scripts/compare_eckardt_diffuser_traverse.py \
  --executable ./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep --allow-cross-family-screening

# Sensitivity to using the traverse-integrated flow column instead
python3 scripts/compare_eckardt_diffuser_traverse.py \
  --executable ./build-eckardt-validation/turbocompgpu_eckardt_diffuser_friction_sweep --allow-cross-family-screening \
  --flow-basis traverse \
  --output-dir benchmarks/results-eckardt-cross-family-screening/radial_o_data_vs_a_model/traverse-flow-sensitivity
```

The script evaluates the existing smooth-wall, Galvas-incidence model at all six exact measurement flows. It writes pointwise predictions and a summary to `benchmarks/results-eckardt-cross-family-screening/radial_o_data_vs_a_model/`. A temporary input uses a dummy total-pressure-ratio field required by the benchmark's map interface; that field is not used to compute a traverse prediction.

## Historical cross-family output (not Impeller A qualification)

Earlier numerical results are retained only for traceability. Because the measurements are for Impeller O and the executable models Impeller A, those residuals do not measure Impeller A model error and must not be used to tune its parameters.
