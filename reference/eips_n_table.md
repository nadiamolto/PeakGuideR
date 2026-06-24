# EIPS elemental isotope-pattern reference table

Reference table used by
[`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md)
to evaluate elemental isotope-pattern support. The EIPS scoring step can
use isotope-pattern evidence from several elements, including Br, Cl, N,
O and S.

Precomputed theoretical lookup table used by
[`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md)
to compare observed elemental isotope-pattern ratios against expected
ratios for different atom counts.

## Usage

``` r
data(eips_n_table)

data(eips_n_table)
```

## Format

A data.frame with precomputed elemental isotope-pattern reference values
used by PeakGuideR.

A data.frame with columns:

- element:

  Element considered for isotope-pattern support: Br, Cl, N, O or S.

- model:

  Model used for the isotope-pattern calculation.

- k:

  Expected isotope peak order or isotope-shift multiplier.

- delta:

  Expected m/z shift from the monoisotopic peak.

- p:

  Isotopic abundance or model probability used in the calculation.

- n:

  Number of atoms of the selected element.

- R_theo:

  Expected theoretical isotope ratio for the given element and atom
  count.

## Source

Internal PeakGuideR reference table.

Internal PeakGuideR theoretical reference table.

## Details

This table contains theoretical values for Br, Cl, N, O and S. It is
used to infer the most compatible number of atoms of a given element
from the observed isotope-pattern ratio.

In this object, the `n` in `eips_n_table` refers to the number of atoms
of the selected element, not specifically to nitrogen. The table
contains theoretical element-count reference values for Br, Cl, N, O and
S.
