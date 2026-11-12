# TurboCompGPU References

This directory contains external scientific references,
experimental datasets, geometry specifications, and
validation materials.

## Organization

- nasa-hecc/: NASA High Efficiency Centrifugal Compressor

The uncopied method and model references below are citation-only. Standards
and textbooks remain under their publishers' terms; this repository does not
redistribute them. See
[`docs/METHOD_AND_MODEL_REFERENCE_AUDIT.md`](../docs/METHOD_AND_MODEL_REFERENCE_AUDIT.md)
for the equation-to-code map and the remaining validation work.

## Data management principles

1. Preserve original experimental data without modification.
2. Record the source and SHA-256 checksum of each original file.
3. Store transformed data separately from original data.
4. Document units, measurement stations, and assumptions.
5. Distinguish experimental observations from model predictions.
6. Do not calibrate model parameters against validation cases.
7. Respect source licensing and redistribution restrictions.

External datasets are not automatically assumed to be
redistributable.

## Compressor methods and model references

### Test method, efficiency, and uncertainty

- **ASME PTC 10-2022, _Axial and Centrifugal Compressors_.** The current
  compressor performance-test code covers test setup, instrumentation,
  evaluation, reporting, thermodynamic performance, and test uncertainty.
  It is the principal method reference for output definitions and test-to-test
  comparisons in this project. The standard is commercially published and is
  cited, not copied. [Official ASME record](https://www.asme.org/codes-standards/find-codes-standards/axial-and-centrifugal-compressors).
- **ISO 5389:2005, _Turbocompressors — Performance test code_.** ISO records
  the 2005 edition as current after confirmation in 2022. Its stated scope is
  performance-test preparation, procedure, evaluation, and assessment for
  turbocompressors. Cited, not copied.
  [Official ISO record](https://www.iso.org/standard/39665.html).
- **ASME PTC 19.1-2018 (R2024), _Test Uncertainty_.** Method reference for
  evaluating uncertainty in measurements, parameters, and methods, and
  propagating it to derived results. Cited, not copied.
  [Official ASME record](https://www.asme.org/codes-standards/find-codes-standards/test-uncertainty).
- Matt Taher and Fred Evans, “Centrifugal Compressor Polytropic
  Performance—Improved Rapid Calculation Results—Cubic Polynomial Methods,”
  *International Journal of Turbomachinery, Propulsion and Power* 6(2), 15
  (2021), [doi:10.3390/ijtpp6020015](https://doi.org/10.3390/ijtpp6020015).
  The article develops a piecewise cubic temperature-entropy path for real-gas
  polytropic performance. It is an open-access article with a
  noncommercial/no-derivatives license; we cite it and do not include its PDF
  or reproduce its tables.

### Meanline and component-model background

- F. J. Wiesner, “A Review of Slip Factors for Centrifugal Impellers,”
  *Journal of Engineering for Power* 89(4), 558–566 (1967). The review surveys
  slip-factor methods and compares them with experimental/potential-flow
  results. [DOI 10.1115/1.3616734](https://doi.org/10.1115/1.3616734). The
  angle convention must be checked when mapping any Wiesner form into this
  code.
- Ronald H. Aungier, _Centrifugal Compressors: A Strategy for Aerodynamic
  Design and Analysis_, ASME Press (2000). Broad meanline design and analysis
  reference for velocity triangles, component models, loss treatments, and
  stage performance. [Official ASME book record](https://www.asme.org/publications-submissions/books/find-book/centrifugal-compressors-strategy-aerodynamic-design-analysis).
- David Japikse, _Centrifugal Compressor Design and Performance_, Concepts
  ETI (1996), ISBN 0-933283-03-2. Design/performance reference for impellers,
  diffusers, and stage behavior. [Catalog record](https://library.kaist.ac.kr/search/ctlgSearch/posesn/view.do?bibctrlno=512091&ty=B).
- S. Larry Dixon and Cesare Hall, _Fluid Mechanics and Thermodynamics of
  Turbomachinery_, 7th ed., Butterworth-Heinemann (2014), ISBN
  978-0-12-415954-9. General turbomachinery thermodynamics and centrifugal
  compressor background. [Official Elsevier record](https://shop.elsevier.com/books/fluid-mechanics-and-thermodynamics-of-turbomachinery/dixon/978-0-12-415954-9).

These textbooks support derivation and implementation review; they are not
validation datasets. Correlation-specific numerical comparisons remain tied to
the papers and machine data listed above.
