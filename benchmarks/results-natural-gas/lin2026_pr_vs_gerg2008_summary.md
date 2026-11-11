# Lin-mixture PR versus GERG-2008 screen

This is a fresh comparison of the stored 20-point Lin et al. temperature and
pressure grid. GERG-2008 values were evaluated with `pyaga8` 0.1.18, which wraps
the NIST AGA 8 GERG-2008 implementation. The mixture uses the eight mole
fractions in Lin et al. Table 1. Pressure is converted from Pa to the binding's
kPa input; density, enthalpy, and heat capacities are converted from molar to
mass-specific units.

Enthalpy is shifted separately for each EOS to that EOS's value at 298.15 K and
101.325 kPa, so the pointwise comparison uses a common state reference. The
reported specific heat ratio is `cp/cv`. GERG's isentropic exponent is a
different quantity and is retained in the pointwise CSV for clarity.

| Property | PR maximum absolute deviation | PR mean absolute deviation | Lin paper maximum |
|---|---:|---:|---:|
| Compressibility factor Z | 2.914% | 2.659% | 1.09% |
| Density | 3.004% | 2.735% | 1.11% |
| Enthalpy relative to common reference | 19.512% | 13.209% | 1.35% |
| Specific heat ratio cp/cv | 2.251% | 1.998% | 1.72% |
| cp (diagnostic; not reported in Lin's table) | 0.723% | 0.373% | — |
| cv (diagnostic; not reported in Lin's table) | 2.906% | 2.055% | — |

The PR screen exceeds all four maximum-deviation values reported by Lin et al.
This does **not** reproduce or contradict their private pointwise comparison:
the paper publishes aggregate errors only, and this implementation uses zero
binary interaction coefficients plus constant ideal-gas component heat
capacities because their parameter values are not published. The result does
show that the current mixture configuration is not qualified against GERG-2008
and must not be presented as a passing natural-gas EOS gate. The discrepancy is
not evidence by itself of a coding defect in the PR cubic; missing/approximate
mixture parameters and caloric assumptions need to be resolved before assigning
cause.

Reproduction:

```sh
g++ -std=c++20 -O2 -I "$THERMOGPU_INCLUDE" -I benchmarks \
  benchmarks/natural_gas_property_grid.cpp \
  "$THERMOGPU_SOURCE/src/thermomodel/peng_robinson_mixture.cpp" \
  -o /tmp/natural-gas-property-grid

/tmp/natural-gas-property-grid \
  > benchmarks/results-natural-gas/lin2026_pr_property_grid.csv

python3 scripts/compare_lin2026_pr_to_gerg.py \
  benchmarks/results-natural-gas/lin2026_pr_property_grid.csv \
  --output benchmarks/results-natural-gas/lin2026_pr_vs_gerg2008.csv
```

The comparison requires `pyaga8` installed in the active Python environment.
The pointwise values are in `lin2026_pr_vs_gerg2008.csv`; the unshifted PR
enthalpy is retained alongside the common-reference enthalpy in the PR grid.

Sources: Lin et al. (2026), Table 1–2; NIST AGA 8 Part 2 / GERG-2008
implementation and documentation: <https://pages.nist.gov/AGA8/> and
<https://github.com/usnistgov/AGA8>.
