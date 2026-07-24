# Union detected adducts per neutral mass id

Combines one or more data.frames (each with at least `neutral_mass_id`
and `adduct` columns) into a single lookup of which adducts were
detected computationally for each `neutral_mass_id`. Used to build the
`detected_adducts_set` consumed by
[`compute_standard_adduct_recovery()`](https://nadiamolto.github.io/PeakGuideR/reference/compute_standard_adduct_recovery.md).

## Usage

``` r
compute_detected_adducts_map(...)
```

## Arguments

- ...:

  Data.frames with `neutral_mass_id` and `adduct` columns, or `NULL`
  (ignored).

## Value

A named list keyed by `as.character(neutral_mass_id)`, each element a
character vector of unique detected adducts.
