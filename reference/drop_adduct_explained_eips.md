# Drop adduct-explained elements from per-feature EIPS evidence

[`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)
/
[`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md)
detect elemental isotope patterns from mass difference and spatial
correlation alone, with no knowledge of the adduct assigned to the same
feature elsewhere in the pipeline. When the detected element is one the
assigned adduct itself introduces (see
[`adduct_introduced_eips_element()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_introduced_eips_element.md)),
the pattern comes from the adduct, not the neutral molecule, and must
not count as EIPS evidence for the neutral-mass candidate.

This removes such elements from `eips_elements` on a per-feature basis.
A feature left with no elements loses its EIPS support
(`has_eips = FALSE`, `eips_score = NA`). Evidence for elements the
adduct does not introduce - any N/O/S evidence, or Cl/Br evidence on a
feature assigned a different adduct - is returned unchanged.

## Usage

``` r
drop_adduct_explained_eips(adduct, has_eips, eips_elements, eips_score)
```

## Arguments

- adduct:

  Character vector of assigned adduct names, one per feature.

- has_eips:

  Logical vector, the feature's raw EIPS flag.

- eips_elements:

  Character vector of `";"`-separated EIPS elements, `NA` where the
  feature has no EIPS evidence.

- eips_score:

  Numeric vector, the feature's raw EIPS score.

## Value

A data.frame with adjusted `has_eips`, `eips_elements` and `eips_score`,
one row per input feature.
