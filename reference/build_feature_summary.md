# Build feature-level summary from PeakGuideR relations

Builds a feature-level summary from a PeakGuideR relation table.

The output contains one row per feature and summarizes whether each peak
is involved in isotope, EIPS or adduct relationships.

## Usage

``` r
build_feature_summary(relation_table, adduct_fam = NULL, pkm = NULL)
```

## Arguments

- relation_table:

  Output from
  [`build_relation_table()`](https://nadiamolto.github.io/PeakGuideR/reference/build_relation_table.md).

- adduct_fam:

  Optional output from
  [`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md).

- pkm:

  Optional peak matrix object with `mass`. If supplied, all features in
  `pkm$mass` are included. If `NULL`, only features present in
  `relation_table` are returned.

## Value

A data.frame with one row per feature.
