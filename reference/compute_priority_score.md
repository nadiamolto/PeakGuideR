# Compute a weighted candidate priority score

Combines per-candidate evidence scores into a single `priority_score` as
a weighted mean. `NA` components are handled in one of two ways
depending on what their absence actually means:

- `standard_adduct_recovery_score` and `eips_evidence_score` are only
  ever evaluable under a structural precondition (standards-library
  membership, and the
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)/C13
  gating in
  [`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md),
  respectively) - `NA` there means the check did not apply, not that it
  was made and failed. These two are excluded from the calculation and
  the remaining weights renormalized to sum to 1, exactly as before.

- `mass_error_score`, `adduct_spatial_score`, `isotope_evidence_score`
  and `family_coherence_score` are evaluable for essentially any
  candidate (a second adduct form, a corroborating family member, or a
  C13 isotope pattern could in principle exist for any organic compound
  with carbon). `NA` there is real negative evidence - it was looked for
  and not found - so it is treated as a literal `0` and its weight stays
  in the average. `mass_error_score` is included in this group only for
  conceptual consistency; in practice it is essentially never `NA`,
  since a candidate only reaches this function after already matching on
  mass.

This means a `broad_db_only` candidate supported only by mass (no
isotope, adduct, family or EIPS evidence) no longer reaches a
`priority_score` near its raw `mass_error_score` - the three zero-filled
evidence components actively pull it down instead of dropping out of the
average.

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
#> [1] 0.6411765 0.3214286
compute_priority_score(scores, weights = default_priority_weights())
#> [1] 0.6411765 0.3214286
```
