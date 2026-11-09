# Project gap-close data and run record

Reviewed: 2026-11-09
Code baseline: `df8cc301fadbeb891a97d87e83c5585563155556`

This note records the current gap review without rewriting the dataset-specific
qualification history in
[`EXTERNAL_COMPRESSOR_QUALIFICATION.md`](EXTERNAL_COMPRESSOR_QUALIFICATION.md).

## Surrogate held-out-speed-line qualification

The reproducible evaluator is
[`qualify_stage_surrogate.py`](../benchmarks/surrogate/qualify_stage_surrogate.py).
It trains a fresh 2→64→64→32→3 network for each leave-one-complete-speed-line-
out fold. Input and log-target scaling are computed from training rows only.
Endpoint folds report speed extrapolation; interior folds report interpolation.
The script records the input CSV SHA256 and writes per-row predictions and a
JSON report.

The qualification has not yet been run on the current boundary grid. It
measures approximation of the Peng–Robinson numerical model, not experimental
compressor accuracy. The evaluator has no acceptance threshold yet; select a
pipeline-use error gate before deployment and preserve edge-fold metrics
separately from interior folds. Dell run instructions are in
[`benchmarks/surrogate/README.md`](../benchmarks/surrogate/README.md).

## Recycle transient test extension

The optional `two_pipe_compressor_coupling` test carries both pipe states over
four quasi-steady compressor steps while changing net export demand from 5.0
to 4.7 to 5.3 to 5.0 kg/s and recycle-valve opening from 0 to 0.25 to 0.5 to
0. It checks each pipe's linepack balance, recycle and compressor-node mass
closure, suction-mixer enthalpy balance, and through-flow accounting. The test
holds speed fixed and applies valve commands directly; it does not model
controller PI memory, compressor shaft dynamics, or actuator lag. The Dell
integration CTest passed on 2026-11-09 in 62.31 seconds. A time-step
refinement study remains before treating the transient as numerically qualified.

The adapter remains isothermal with fixed segment temperature and
compressibility. No measured station transient has been reproduced.

## Public physical guide-vane and surge source

Nan Chen's 2017 TUM dissertation, *Experimental and Numerical Investigation
of a Centrifugal Compressor Stage with Variable Inlet Guide Vanes*, is
available from the university repository:
<https://mediatum.ub.tum.de/doc/1336387/1336387.pdf>.

Its experimental rotating campaign covers three IGV designs at −20°, 0°,
20°, 40°, and 60°, with at least five steady points per measured speed line.
The exact normalized nominal flow coefficients from Table 3.3 and quantitative
comparisons stated in the Figure 3.17 discussion are transcribed in
[`tum_chen_2017_vigv_reported_results.csv`](../data/validation/tum_chen_2017_vigv_reported_results.csv).
The source reports normalized nominal flow coefficients 1.152, 1.000, 0.944,
0.751, and 0.616 for those respective angles. It says B/C pressure ratio is
at least 2% above A at 0° design flow; B/C are nearly equal and about 3% above
A at 20° nominal flow; Type C reaches a gap up to 6% over A at 60° nominal
flow. It reports up to 2% stage-efficiency improvement and about 5% head
coefficient improvement for B/C at 0° design flow. These comparisons retain
the source's qualifiers; they are not treated as exact machine-map points.

The dissertation supplies plotted speed lines rather than a machine-readable
coordinate table. The blank import schema is
[`tum_chen_2017_vigv_template.csv`](../data/validation/tum_chen_2017_vigv_template.csv);
actual per-speed curve coordinates and surge-line digitization still need to
be extracted with uncertainty. The adopted surge procedure throttles in 0.1%
steps every five seconds until dynamic pressure exceeds 100 mbar peak-to-peak.
At 20° the source reports rotating stall before this threshold. The three vane
designs have minimum surge flow within 1% at each angle.

The TUM compressor is different hardware. Its results can test the map-import
and interpolation path, but cannot calibrate the project's Eckardt/Sandia
cases or provide their anti-surge boundary. Machine-specific guide-vane and
surge data are still needed. Confirm reuse terms before committing digitized
curve coordinates.

## Machine-map software boundary

The map-import pass adds complete-grid CSV loading and bounded interpolation
for surge flow, guide-vane effects, recycle-valve flow, and cooler return
temperature. The interfaces preserve nodewise conservative uncertainty bounds;
the surge helper can report a margin against an uncertainty-adjusted boundary.
The recycle return adapter now passes current recycle flow to three-argument
cooler models while preserving the previous two-argument callback interface.
The `calibrated_map_interpolation` unit test covers interpolation, uncertainty,
compatibility, and rejection of out-of-domain queries. It passed on Dell 7710
on 2026-11-09: 1/1 CTest passed in 0.01 seconds. This is a software-interface
check, not physical validation of any compressor map.

These are data interfaces, not populated machine maps. All four map templates
remain header-only, and anti-surge must stay unavailable until a selected-unit
surge line and its event definition are sourced and independently qualified.
See [`CALIBRATED_MAPS.md`](CALIBRATED_MAPS.md) for column names, domains, and
map-specific acquisition requirements.

## Sandia effective-parameter exploratory result

In the Dell experiment, the first bounded fit used a lower bound of 0.60 and
landed on that bound. After widening the search interval to 0.30–0.60, the
fitted effective parameter was 0.367983. On the five points at or below
39,000 rpm, training MAPE was 0.777%. On the seven higher-speed points,
holdout MAPE was 3.142% and maximum absolute error was 7.324%. All seven
higher-speed errors were negative; the largest magnitude was at 56,000 rpm.

Because the high-speed errors were inspected before expanding the bounds,
this is exploratory transfer evidence, not a pristine blind holdout. The
fitted scalar is an empirical head correction, not measured blade slip. The
dataset is a Sandia sCO2 impeller-only screen and does not qualify full-stage
CO2 behavior.
