# CIR score \#' Validate C13 M+1 candidates using carbon isotope ratios

Filters C13 M+1 candidates from
[`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)
and compares the observed M+1/M+0 intensity ratio against the
theoretical carbon isotope ratio expected for the monoisotopic m/z.

## Usage

``` r
cir_score(
  result,
  pkm,
  min_score_final = 0.6,
  cir_rel_tol = 0.3,
  ratio_method = c("sum", "mean", "median"),
  min_quantile = 0.01,
  mask_strategy = c("M0_only", "AND")
)
```

## Arguments

- result:

  Isotope detection results
  ([`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)
  output)

- pkm:

  Peak matrix list (`mass`, `intensity`)

- min_score_final:

  Minimum `score_final` threshold in the range 0 to 1. Default: 0.6

- cir_rel_tol:

  CIR relative tolerance in the range 0 to 1. Default: 0.3

- ratio_method:

  `"sum"` (default), `"median"`, or `"mean"` for I_M1/I_M0

- min_quantile:

  Quantile mask for M0 intensities. Default 0.01

- mask_strategy:

  `"M0_only"` (recommended) or `"AND"`. Default `"M0_only"`.

## Value

A data.frame with one row per C13 M+1 candidate and the following
additional columns: `R_obs`, `R_theo`, `cir_rel_err`,
`is_valid_c13_raw`, `cir_score`, `cir_class`, `is_chained_c13`,
`is_valid_c13`, `has_C13_M2`, `idx_C13_M2`, `mz_C13_M2` and
`score_C13_M2`.

## Details

The function is intended to add isotope-ratio evidence to
morphology-based C13 candidates. It does not assign definitive isotope
annotations by itself.

Filters C13_M1 candidates and validates against theoretical ratios.
