# Recover compound_db candidates that top_n truncated away from a standard-linked neutral mass

`top_n` in
[`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md)
is a blind cut by ascending `candidate_ppm_error`: it has no visibility
into `standards_db`, so a candidate that would otherwise resolve to
`source == "both"` in
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)
can be truncated away simply because a common formula/isobar crowds the
same mass window with more than `top_n` `compound_db` entries within
tolerance - the standards-library match is exactly the strongest
evidence the candidate could have, and it never gets the chance to
compete.

This function undoes that for the (typically very small) subset of
`neutral_mass_id` values where `standard_candidates` found at least one
hit: it re-runs
[`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md)
for those masses only, with `top_n = NULL` (no cut), and unions the
recovered rows into `compound_candidates`, deduplicating on
`neutral_mass_id` + `compound_identity_id` so a candidate that already
survived the original cut is not duplicated. `neutral_mass_id` values
with no standard-library hit are untouched - `top_n` still applies to
the overwhelming majority of masses exactly as before, so this has no
effect on runtime or output size outside the small subset where a
standard is actually in play.

## Usage

``` r
expand_compound_candidates_for_standards(
  compound_candidates,
  standard_candidates,
  neutral_for_matching,
  compound_db,
  ppm_tol
)
```

## Arguments

- compound_candidates:

  Output of
  [`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md),
  possibly truncated by `top_n`.

- standard_candidates:

  Output of
  [`match_standards_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_standards_by_mass.md),
  or `NULL`/empty if standard-adduct support is not being used.

- neutral_for_matching:

  The same data.frame (`neutral_mass_id`, `neutral_mass_consensus`)
  passed to
  [`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md)
  to produce `compound_candidates`.

- compound_db:

  The same `compound_db` passed to
  [`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md).

- ppm_tol:

  The same `ppm_tol` passed to
  [`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md).

## Value

`compound_candidates`, with any `compound_db` candidates that `top_n`
had excluded for a standard-linked `neutral_mass_id` added back in.
