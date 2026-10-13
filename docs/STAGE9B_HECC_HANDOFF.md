# TurboCompGPU — NASA HECC Validation Handoff

Date: 2026-10-13

## Current status

Stage 9A: Complete.
- Frictional impeller and vaneless diffuser integration.
- IdealGas and PengRobinson verification.
- Commit: 013ddee.

Stage 9B.1: References acquired.
- NASA HECC publications and archive documentation.
- Commit: 1a678b4.

Stage 9B.2: Complete.
- Variable-width frictional diffuser.
- Shared RK4 integrator.
- Analytical and finite-difference EOS derivatives.
- Independent frictionless algebraic verification.
- Fourth-order convergence verified for both EOS.
- Commit: a664ffc.
- Regression baseline: 69 enabled tests passing,
  one test intentionally disabled.

Stage 9B.3: Geometry infrastructure.
- PCHIP interpolation for radial hub/shroud profiles.
- Analytic first derivatives.
- Passage width and width derivative.
- Input validation and positive-width checks.
- Manufactured-solution and solver integration tests.
- No NASA geometry has been reconstructed yet.

## NASA archive access blocker

Published URL:
https://storage.googleapis.com/hecc-data/NASA-HECC-Data-Archive.zip

Observed on 2026-10-13:
HTTP 403, Google Cloud Storage AccessDenied.

Response says anonymous caller lacks storage.objects.get
permission, or the object may not exist.

The archive has NOT been downloaded.

The original NASA hot hub/shroud coordinates and experimental
vaneless diffuser CSV files are NOT available locally.

Do not represent synthetic geometry as NASA geometry.

## Geometry reference dimensions

Impeller trailing edge:
- Radius: approximately 0.21542 m
- Blade height: approximately 0.01524 m

Diffuser rating plane:
- Radius: approximately 0.23134 m
- Passage width: approximately 0.01222 m

Stage exit:
- Approximate midspan radius: 0.30545 m
- Passage width: approximately 0.00795 m

These are published reference dimensions, not a complete
reconstruction of the actual hot wall geometry.

## Important modeling assumptions

Current variable-width diffuser uses:

    A(r) = 2*pi*r*b(r)

and:

    d(ln A)/dr = 1/r + b'(r)/b(r)

The model assumes a radial annular passage.
For strongly curved walls, assess the difference between
axial wall separation and flow-normal passage width.

PCHIP is applied independently to hub and shroud.
The resulting width is not necessarily shape-preserving,
so positivity must be checked independently.

Coordinate units and conventions must be established
from NASA source data before ingestion.

## Next steps

1. Obtain original NASA HECC archive or coordinate files.
2. Confirm source units and axial/radial conventions.
3. Implement a strict coordinate-file reader with provenance.
4. Load actual hot hub/shroud coordinates.
5. Verify stations and published dimensions.
6. Evaluate radial-area approximation for the actual walls.
7. Connect measured geometry to the verified diffuser solver.
8. Compare predictions against NASA experimental data.

## Resume instructions

Repository: rdm375/TurboCompGPU
Branch: main

    cd ~/Sync/TurboCompGPU
    git status --short
    git log -5 --oneline
    cat docs/STAGE9B_HECC_HANDOFF.md
    cmake --build build-validation -j2
    ctest --test-dir build-validation --output-on-failure

Do not claim NASA experimental validation until the original
geometry and measurements have been obtained and analyzed.
