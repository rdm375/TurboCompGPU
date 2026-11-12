# GPU friction closure

The batched CUDA backend evaluates the Peng–Robinson mixture state and the
friction closure on the GPU. In the vaneless diffuser, each Runge–Kutta stage
evaluates the local PR density and velocity, sets the annular hydraulic
diameter to `Dh = 2 b`, computes Reynolds number, and calls CUDA Colebrook with
the operating point's relative roughness. With 512 radial steps, this is 2,048
local Colebrook evaluations per operating point. CPU parity uses the same
stage-local equations and the CPU Colebrook implementation.

At each diffuser RK stage the default path evaluates natural-gas dynamic
viscosity with the Lohrenz–Bray–Clark (LBC) method specified in ISO 20765-5:2022
§5.1. It uses the fixed mixture composition, local temperature and pressure,
and local PR mass density converted to molar density. Relative roughness stays
constant. Thus density, velocity, viscosity, Reynolds number, and Darcy factor
are local, while wall roughness is a fixed geometry input. The API retains a
constant-viscosity fallback for comparisons.

The ISO method expects density from ISO 20765-2 / GERG-2008, so this
implementation uses PR density as a screening substitution. It must not be
described as ISO-conforming, and the standard's approximately 4% viscosity
uncertainty estimate does not directly apply to this PR-density substitution.
The implementation currently enforces the standard's stated gas metering
range (253.15–313.15 K and 0.1–10 MPa) and supports components for which the
ISO 20765-2 critical molar density table is encoded.

Impeller passage friction remains a lumped Darcy–Weisbach closure. It uses
inlet relative velocity, inlet density, the same ISO LBC viscosity evaluated
at the inlet state (or the configured fallback), passage length,
and hydraulic diameter. The input API does not yet contain blade-to-blade or
meridional impeller stations from which a local passage state and local friction
factor can be integrated. Its output Darcy factor describes that impeller
screening estimate; it is not the diffuser's mean factor.

Impeller passage friction remains a lumped Darcy–Weisbach closure using inlet
relative speed and inlet viscosity. A distributed impeller friction model
requires a geometry/state path through the impeller; it should be added when
those station definitions are part of the input model, not inferred from a
single inlet/outlet pair.

References: ISO 20765-5:2022, §5.1 (LBC viscosity method), and ISO 20765-2:2015,
Annex B (component critical properties). The standard notes density sensitivity
and recommends GERG-2008 density for the LBC calculation.
