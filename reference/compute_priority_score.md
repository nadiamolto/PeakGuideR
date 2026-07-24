# Compute a weighted candidate priority score

Combines per-candidate evidence scores into a single `priority_score` as
a weighted mean, excluding `NA` components from the calculation and
renormalizing the remaining weights to sum to 1. `NA` components are
never treated as `0`.

## Usage

``` r
compute_priority_score(scores, weights = default_priority_weights())
```

## Arguments

- scores:

  A data.frame with one row per candidate and columns named after
  `weights` (missing columns are treated as all-`NA`), or a single named
  numeric vector for one candidate.

- weights:

  Named numeric vector of weights, see
  [`default_priority_weights()`](https://nadiamolto.github.io/PeakGuideR/reference/default_priority_weights.md).

## Value

A numeric vector (or scalar) of priority scores in the range 0 to 1, or
`NA` where every component is `NA`.

## Examples

``` r
scores <- data.frame(
  mass_error_score = c(0.95, 0.80),
  adduct_spatial_score = c(0.70, NA),
  isotope_evidence_score = c(NA, NA),
  eips_evidence_score = c(NA, NA),
  standard_adduct_recovery_score = c(1, NA),
  family_coherence_score = c(0.65, 0.65)
)

compute_priority_score(scores)
#> [1] 0.8384615 0.7500000
compute_priority_score(scores, weights = default_priority_weights())
#> [1] 0.8384615 0.7500000
```
