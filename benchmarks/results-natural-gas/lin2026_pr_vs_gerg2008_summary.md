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
| Enthalpy relative to common reference | 19.512% | 13.209% | Not directly comparable* |
| Specific heat ratio cp/cv | 2.251% | 1.998% | 1.72% |
| cp (diagnostic; not reported in Lin's table) | 0.723% | 0.373% | — |
| cv (diagnostic; not reported in Lin's table) | 2.906% | 2.055% | — |

The PR screen exceeds the Lin paper's reported maximum deviations for Z,
density, and cp/cv. The enthalpy percentage is not a like-for-like score: the
paper does not state its enthalpy reference datum, while this comparison shifts
each EOS to a common 298.15 K, 101.325 kPa state. Since enthalpy itself has an
arbitrary additive reference, the 19.512% figure must not be compared directly
with the paper's 1.35% limit. The matched-reference values still permit
reference-invariant enthalpy-difference checks. Across each isotherm, the
6.5-to-8.0 MPa enthalpy-change discrepancy is 1.00–1.06 kJ/kg, or 5.9–7.1%
of the GERG enthalpy change; this is a property-grid diagnostic, not a
compressor-work prediction.

This does **not** reproduce or contradict the paper's private pointwise
comparison: Lin et al. publish aggregate errors only, and this implementation
uses zero binary interaction coefficients plus constant ideal-gas component
heat capacities because their parameter values are not published. The three
directly comparable failures mean this PR mixture is not qualified against
GERG-2008 and must not be presented as a passing natural-gas EOS gate. The
discrepancy is not evidence by itself of a coding defect in the PR cubic;
missing/approximate mixture parameters and caloric assumptions need to be
resolved before assigning cause.

\*The paper reports a 1.35% maximum deviation for specific enthalpy, but does
not document the enthalpy reference convention needed to compare that
percentage with the common-reference result here.

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
