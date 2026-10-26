# Append the mass defect flag to a candidate table

Adds a `mass_defect_ok` column at the end of a table with a
`neutral_mass_consensus` column (`neutral_mass_candidates` or
`candidate_annotations`), using
[`mass_defect_ok()`](https://nadiamolto.github.io/PeakGuideR/reference/mass_defect_ok.md).
The remaining columns and rows are left untouched. If the table has no
`neutral_mass_consensus` column, it is returned unchanged.

## Usage

``` r
add_mass_defect_flag(x)
```

## Arguments

- x:

  A data.frame with a `neutral_mass_consensus` column.

## Value

`x` with the extra column `mass_defect_ok`.
