# Calibrated guide-vane, surge, and recycle maps

Reviewed: 2026-11-09

This document defines the tabular calibration interface added for the
compressor and recycle maps. It does not provide selected-machine values.
Every template under `data/validation/*_map_template.csv` is intentionally
header-only. Do not enable machine control from a template, a synthetic test
fixture, or the cross-machine TUM summary.

## Interpolation contract

`RectilinearCalibrationTable<N, M>` loads a complete rectilinear grid from CSV
and multilinearly interpolates `M` outputs over `N` numeric axes. It validates
finite values, monotone axes, unique grid nodes, and a complete tensor grid.
Queries on the boundary are allowed; queries outside any calibrated axis throw
`std::out_of_range`. It never extrapolates or silently clamps a query.

Example map headers and loaders:

```cpp
using SurgeTable = turbocompgpu::SurgeLineTable;
auto surge_table = SurgeTable::load_csv(
    "surge.csv",
    {"corrected_speed_rpm", "guide_vane_angle_deg"},
    {"surge_corrected_flow_kg_s", "surge_flow_uncertainty_bound_kg_s"});
turbocompgpu::TabulatedSurgeLineMap surge(std::move(surge_table));
const double margin = surge.margin_fraction(
    corrected_speed_rpm, guide_vane_angle_deg,
    corrected_pipeline_flow_kg_s);
const double conservative_margin = surge.conservative_margin_fraction(
    corrected_speed_rpm, guide_vane_angle_deg,
    corrected_pipeline_flow_kg_s, 1.0);
```

The returned margin is `(corrected flow - surge flow) / surge flow`; positive
is on the higher-flow side, zero is on the mapped line, and negative is below
the line. Speed and flow corrections must use the same reference conditions as
the data source. The caller supplies this qualified margin to the controller;
the controller does not infer the correction or map.

The map linearly interpolates its second output. Supply a conservative
absolute-error bound at each node, rather than a one-sigma standard deviation,
unless the interpolation is changed to propagate covariance. The conservative
helper moves the boundary upward by the requested multiple of that bound
before calculating margin. It does not replace a qualification policy or
account for correlated calibration uncertainty.

## Guide-vane effects table

`guide_vane_effects_map_template.csv` uses the axes, in order,
`guide_vane_angle_deg`, `speed_rpm`, `compressor_flow_kg_s`; the four outputs
are `inlet_prewhirl_m_s`, `entropy_generation_j_kg_k`,
`prewhirl_uncertainty_bound_m_s`, and
`entropy_uncertainty_bound_j_kg_k`. The uncertainty outputs are linearly
interpolated conservative bounds. Load it as
`GuideVaneEffectsTable`, then call
`TabulatedGuideVaneCalibration(table).bind(domain)` to obtain the existing
bounded `GuideVaneMap` consumed by `CompressorComponent`.

This is a fitted effects map at the stated conditions. It must not interpret a
vane command angle as delivered swirl without measured or otherwise qualified
turning data. If inlet state changes materially, calibrate a corrected map or
reject queries outside the validated conditions; this first table form has no
inlet-state axes.

## Recycle valve and cooler-return tables

`recycle_valve_map_template.csv` axes, in order, are `opening`,
`suction_pressure_pa`, `discharge_pressure_pa`, and
`discharge_temperature_k`; outputs are measured `flow_kg_s` and its
conservative uncertainty bound `flow_uncertainty_bound_kg_s`. Bind a loaded
`RecycleValveFlowTable` through `TabulatedRecycleValveCalibration` and set the
corresponding explicit `RecycleValveDomain`, including a measured maximum flow.
The solver still enforces the selected EOS state and pressure direction.

`recycle_cooler_return_map_template.csv` axes, in order, are
`discharge_pressure_pa`, `discharge_temperature_k`, `suction_pressure_pa`, and
`recycle_flow_kg_s`; outputs are `return_temperature_k` and conservative
uncertainty bound `temperature_uncertainty_bound_k`. Bind a loaded
`RecycleCoolerReturnTable` to the chosen EOS model. Its three-argument return
callback uses the current valve flow, reconstructs the EOS state at suction
pressure, and lets the existing energy balance report cooler duty. Historical
two-argument return callbacks remain supported, but they cannot represent a
flow-dependent cooler response. Query
`return_temperature_uncertainty_k(...)` separately when reporting the
calibration uncertainty.

The return table represents the combined measured recycle return path and
cooler at the tabulated boundary states. If pressure loss, cooler bypass, or
valve hysteresis is important, model and calibrate those separately rather
than hiding them in a single temperature curve.

## Calibration and qualification required

For each map, preserve source, machine/serial, gas composition, reference
conditions, measurement uncertainty, calibration date, and units. The current
interpolators treat uncertainty outputs as absolute bounds and interpolate
them linearly; do not enter one-sigma values without a suitable propagation
method. Record
whether the axes are actual or corrected quantities. Keep train/calibration
points separate from withheld validation points. Inspect the map edges and
uncertainty before binding it to a component.

- Surge: measure the stability boundary over corrected speed and every
  supported vane position using a documented surge/stall event definition.
  Do not substitute the last stable solver point for a measured surge point.
- Guide vane: measure delivered inlet swirl/prewhirl and total-pressure loss
  (or a justified entropy equivalent) over angle, speed, and flow. Check angle
  sign and actual vane travel.
- Recycle valve: record flow over opening and pressure/temperature conditions
  with the actual gas or a documented conversion. Include shutoff and maximum
  travel behavior.
- Cooler/return: record return temperature over flow and inlet/outlet
  conditions. Check outlet state and duty with the same EOS used in the
  compressor calculation.

The TUM summary data in `tum_chen_2017_vigv_reported_results.csv` provide
cross-machine trends only. They do not fill these map tables. Until selected-
unit maps are loaded and independently checked, keep anti-surge disabled and
treat guide-vane/recycle tables as unavailable for operational use.
