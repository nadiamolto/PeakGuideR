# Build neutral-mass candidate table

Builds a neutral-mass candidate table from adduct families.

This function summarizes neutral masses inferred from adduct families,
adds isotope and EIPS support from `feature_summary`, matches the
inferred neutral masses against a compound mass database, and optionally
checks whether each candidate compound and inferred adduct are supported
by the standard-adduct library.

The output contains one row per inferred neutral mass and compound
candidate. Candidate compounds are putative mass matches, not definitive
identifications.

## Usage

``` r
build_neutral_mass_candidates(
  adduct_fam,
  feature_summary,
  compound_db = NULL,
  standards_db = NULL,
  ion_mode = c("pos", "neg"),
  matrix = NULL,
  ppm_tol = 5,
  neutral_cluster_ppm = 5,
  top_n = 10L,
  quiet = FALSE
)
```

## Arguments

- adduct_fam:

  Output from
  [`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md).

- feature_summary:

  Output from
  [`build_feature_summary()`](https://nadiamolto.github.io/PeakGuideR/reference/build_feature_summary.md).

- compound_db:

  Compound mass database. If `NULL`, the included non-commercial
  compound mass database is loaded.

- standards_db:

  Standard adduct library. If `NULL` and `matrix = "HCCA"`, the included
  non-commercial standard adduct library is loaded.

- ion_mode:

  `"pos"` or `"neg"`.

- matrix:

  Matrix name. Standard-adduct support is currently applied only when
  `matrix = "HCCA"`. Use `NULL` to skip standard-adduct support.

- ppm_tol:

  PPM tolerance for compound mass matching.

- neutral_cluster_ppm:

  PPM tolerance used to group similar neutral masses inferred from
  different adduct families.

- top_n:

  Maximum number of compound candidates kept per neutral mass.

- quiet:

  Logical. If `FALSE`, database loading functions may print notices.

## Value

A data.frame with one row per neutral mass and compound candidate.
