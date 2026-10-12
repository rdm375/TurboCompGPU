# Stage 7J — Outlet Geometry Sensitivity Validation

## Scope

Numerical investigation of compressor characteristic sensitivity
to impeller outlet radius and width.

The study uses the existing Peng–Robinson methane model,
incidence and diffusion losses, and independently evaluated
provisional passage-turning choking capacity.

Characteristic maxima are not validated surge boundaries.
The provisional capacity is not a fully coupled choking solution.

## Stage 7J — Independent Geometry Sensitivity

45 speed/geometry combinations.

Baseline speeds:

- 9000 RPM
- 10500 RPM
- 12000 RPM
- 13500 RPM
- 15000 RPM

Outlet radius and width were perturbed independently.

The provisional choking capacity remained invariant under
these outlet-geometry perturbations in the implemented model.

## Stage 7J.1 — Outlet Width Transition

At 15000 RPM:

Critical outlet width:

    26.404273808 mm

Relative increase:

    5.617095232%

Transition determined from the zero of the upper-boundary
pressure-ratio derivative.

## Stage 7J.2 — Constant-Area Geometry Sensitivity

Transformation:

    r2(s) = s * r2_baseline
    b2(s) = b2_baseline / s

Effective outlet area is invariant.

The local elasticity identity is:

    E_constant_area = E_radius - E_width

The numerical results agree with this identity within the
finite-difference truncation errors of the sensitivity estimates.

## Stage 7J.3 — Constant-Area Boundary Transition

At 15000 RPM:

    critical_scale  = 0.963718071207
    critical_radius = 192.743614241 mm
    critical_width  = 25.941196650 mm

The root converged under refinement of the mass-flow
derivative step.

The finest-three-step scale spread was:

    1.013278961e-7

The maximum relative outlet-area error was:

    2.220446049e-16

## Constant-Area Peak-Flow Elasticities

Richardson-extrapolated estimates:

| RPM | Elasticity |
|---:|---:|
| 9000 | -4.1547722757 |
| 10500 | -2.9425492371 |
| 12000 | -2.3248887173 |
| 13500 | -1.9230788201 |
| 15000 | -1.6226751545 |

The 15000 RPM symmetric finite difference is unavailable
at a 4% perturbation because the lower geometry crosses
the provisional capacity boundary.

Smaller perturbations remain on the interior-peak branch.

The 9000 RPM convergence sequence exhibits some departure
from ideal second-order behavior at the finest step.

## Limitations

- Provisional choking capacity is computed independently.
- No coupled characteristic/choking closure is claimed.
- No experimentally validated surge boundary is claimed.
- Finite-difference convergence is not a physical uncertainty bound.
- Similarity-coordinate superiority remains to be demonstrated.

## Next Direction

Investigate compressor-map representations using physically
motivated transformed coordinates, including:

- Effective outlet flow area.
- Impeller tip speed.
- Flow coefficient.
- Head coefficient.
- Surge-to-choke normalized coordinates where appropriate
  and supported by validated boundary definitions.

Evaluate these transformations for adaptive surrogate sampling,
approximation error, and GPU execution efficiency.
