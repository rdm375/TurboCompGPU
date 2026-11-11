# Eckardt B speed-line holdout

**Scope:** this is a same-family, structured speed-line holdout. The
effective slip/head parameter is fitted on the listed training lines.
It is empirical and must not be interpreted as a physical slip
measurement or a cross-family qualification.

- Training speed lines: 10000, 12000, 14000 rpm (13 points)
- Held-out speed lines: 16000 rpm (4 points)
- Fitted effective parameter: 0.74572749
- Training MAPE: 0.921%
- Held-out MAPE: 3.785%
- Held-out maximum absolute error: 5.151%
- Predeclared held-out maximum-error gate: 7.000%
- Result: **PASS**

The 7% default is a provisional project gate informed by the maximum
Eckardt B pressure-ratio deviation reported for a separate published
mean-line model. It is not an experimental uncertainty bound. Change
the gate only before fitting, and record the reason.

Pointwise results: `/tmp/turbocompgpu-remote-physics/eckardt-b-holdout.csv`.
