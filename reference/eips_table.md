# EIPS theoretical lookup table

Precomputed expected isotope presence scores (EIPS) for different
isotope types and atom counts.

Database-derived lookup table used by
[`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md)
to evaluate optional formula support for elemental isotope-pattern
evidence.

## Usage

``` r
data(eips_table)

data(eips_table)
```

## Format

A data.frame with columns used to store theoretical or empirical isotope
presence score expectations.

A data.frame with columns:

- formula:

  Molecular formula.

- mz_mono:

  Monoisotopic mass or monoisotopic m/z reference value.

- element:

  Element considered for isotope-pattern support: Br, Cl, N, O or S.

- model:

  Model used for the isotope-pattern calculation.

- k:

  Expected isotope peak order or isotope-shift multiplier.

- delta:

  Expected m/z shift from the monoisotopic peak.

- n_el:

  Number of atoms of the selected element in the formula.

- R_theo:

  Expected theoretical isotope ratio for the formula-level entry.

## Source

Internal PeakGuideR reference table.

Internal PeakGuideR reference table derived from the compound mass
database distributed with PeakGuideR.

## Details

This table links candidate formulas and monoisotopic masses to the
number of atoms of selected non-carbon elements. It is used to check
whether an inferred neutral mass and estimated element count are
compatible with database-derived formula candidates.

This object is used by
[`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md)
only when formula support is evaluated. It is joined by `element` and
filtered by neutral-mass agreement and agreement between the inferred
element count and the formula-derived element count.
