# Eckardt A efficiency cross-check

## Data source and limits

`data/validation/eckardt_impeller_a_efficiency_digitized.csv` is an approximate
digitization of the measured symbols in panel (d) of Figure 6 in Meroni et al.,
“Design of centrifugal compressors for heat pump systems,” *Applied Energy*
232 (2018). The peer-reviewed postprint is available from the
[DTU repository](https://backend.orbit.dtu.dk/ws/files/154611324/Postprint_Design_of_centrifugal_compressors_for_heat_pump_systems.pdf).
No raw table or pointwise uncertainty was found. Estimated digitization
uncertainty is about 0.05 kg/s in flow and 0.01 absolute in efficiency. These
are screening values, not precision measurements. Flow points are digitized
independently of the Bayomi pressure-ratio transcription and are not a
point-for-point paired dataset.

## Method

The benchmark computes total-to-total efficiency from its calculated total
states:

`eta_tt = (T03s - T01) / (T03 - T01)`, where `T03s` comes from the pressure
ratio `p03/p01` using the constant-cp ideal-gas isentropic relation. The
impeller model has no aerodynamic loss term, so diffuser friction is the only
loss mechanism in this efficiency calculation.

## Result

The later smooth-wall inlet-Reynolds model with Galvas incidence completed all
21 digitized points. Its mean signed and mean absolute efficiency errors were
both **+4.47 percentage points**; maximum absolute error was **23.91 points**.
The no-incidence variant gave **+5.18 points** mean signed/absolute error and
**23.98 points** maximum absolute error. Thus the incidence correction makes
only a small change to aggregate efficiency error; high-flow residuals remain
large. Earlier lossless/impeller-loss-free screening numbers are superseded by
this model variant.

The pressure-ratio sweeps show that diffuser wall friction and impeller-loss
terms both affect agreement, but a good pressure-ratio fit does not guarantee
good efficiency. The efficiency residual remains large at some high-flow
points. No friction factor was selected to fit either dataset.

This is a screening comparison, not qualification. The max error is sensitive
to digitized points at high flow, and the simplified smooth-wall correlation is
held constant along the diffuser. Meroni et al. report efficiency discrepancies
up to 6.7 percentage points at high flow for their own mean-line model; that
result is context, not a direct apples-to-apples target for this model.

## Paired pressure/work audit

The pressure-ratio transcription and efficiency digitization are not paired
measurements. To see how much of the efficiency/work discrepancy is associated
with pressure ratio, `scripts/audit_eckardt_efficiency_energy.py` linearly
interpolates the measured pressure-ratio speedline at each digitized efficiency
flow. It does not extrapolate. In the current data, 17 of 21 efficiency samples
are inside their corresponding pressure speedline; four are explicitly marked
`outside_measured_speedline` and receive no derived temperature or work.

For paired points, the script computes ideal-gas `T03s` from the interpolated
pressure ratio and obtains an implied `T03` and specific work from the digitized
efficiency. It computes the same quantities from the model's predicted
pressure ratio and efficiency. The work residual is then split exactly into
two ordered contributions:

`pressure-ratio contribution = (dh_is,pred - dh_is,meas) / eta_meas`

`efficiency contribution = dh_is,pred * (1 / eta_pred - 1 / eta_meas)`

These sum to predicted-minus-implied specific work. This is an accounting
decomposition, not proof of physical root cause; its component values depend
on the stated order of the split. It uses constant `cp=1004.5 J/(kg K)`,
`R=287.05 J/(kg K)`, and `T01=288.15 K`, matching the ideal-gas screening
assumptions. Approximate efficiency digitization and the non-paired source data
remain the main limitations.

On the 17 paired points, the smooth-wall plus Galvas-incidence rows have
interpolated pressure-ratio MAPE **2.539%** and mean work residual
**−6.776 kJ/kg**. The ordered decomposition attributes **−2.942 kJ/kg** to
pressure-ratio difference and **−3.834 kJ/kg** to the efficiency difference.
The smooth-wall no-incidence comparison gives **2.257%** pressure-ratio MAPE
and **−6.742 kJ/kg** mean work residual, split into **−2.609** and **−4.133
kJ/kg**, respectively. Negative residual means the model predicts less work
than implied by the digitized efficiency and interpolated pressure ratio.

These figures indicate that pressure-ratio underprediction contributes to the
work gap, while the predicted efficiency also remains too high. They do not
identify which omitted loss mechanism causes the efficiency error. Run the
audit after generating the incidence/choke efficiency comparison:

```sh
python3 scripts/audit_eckardt_efficiency_energy.py
```

The default per-point output is
`benchmarks/results-eckardt-impeller-a/paired_energy_audit.csv`. It records the
pressure interpolation bracket, paired/unpaired status, implied and predicted
outlet temperature/work, residual components, and model variants. Alternate
input/output paths and gas assumptions can be supplied through `--help`.
