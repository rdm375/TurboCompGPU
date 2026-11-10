# Eckardt A distributed skin-friction screen

This diagnostic compares the current Jansen mean-state friction closure with
an integral using the same smooth-wall Fanning factor at local Reynolds
numbers along a one-dimensional passage.

The screen integrates `dh_sf = 2 Cf(Re) W^2 ds / Dh` with composite Simpson
quadrature over 64 intervals. Relative speed, density, and temperature vary
linearly from the mean inlet state to the exit state. The Galvas equivalent
hydraulic diameter and Jansen path length are held fixed. This isolates the
mean-Re approximation; it does not validate the assumed state path or passage
geometry, and it does not replace the baseline closure.

Compare the `distributed_skin_friction_screen` row with
`all_losses_enthalpy_equivalent` in `loss_term_ablation.csv`. The Reynolds
number and Fanning-factor columns remain the baseline mean-state diagnostics;
the scenario's integrated loss is in `impeller_skin_friction_loss_j_kg`.
Treat any improvement on Eckardt A as a sensitivity result, not independent
qualification.

## Geometry-continuity diagnostic

`geometry_continuity_skin_friction_screen` interpolates effective inlet/exit
flow area and rectangular endpoint hydraulic diameters. It obtains local
meridional speed from continuity, then combines that with linearly interpolated
radius and absolute swirl to form relative speed. Density and temperature are
also interpolated. The path length and smooth-wall friction law stay fixed.

This is a kinematic sensitivity screen, not an energy-closed throughflow
solution. Compare it with `all_losses_enthalpy_equivalent` and
`distributed_skin_friction_screen`; do not treat an improved fit as validation
of the assumed internal path.
