# Default weights for the candidate priority score

Returns the default weights used by
[`compute_priority_score()`](https://nadiamolto.github.io/PeakGuideR/reference/compute_priority_score.md)
/
[`build_candidate_annotations()`](https://nadiamolto.github.io/PeakGuideR/reference/build_candidate_annotations.md)
to combine per-candidate evidence scores into a single `priority_score`.

## Usage

``` r
default_priority_weights()
```

## Value

A named numeric vector of weights summing to 1: `mass_error_score`
(0.20), `adduct_spatial_score` (0.20), `isotope_evidence_score` (0.20),
`eips_evidence_score` (0.15), `standard_adduct_recovery_score` (0.15),
`family_coherence_score` (0.10).

## Examples

``` r
default_priority_weights()
#>               mass_error_score           adduct_spatial_score 
#>                           0.20                           0.20 
#>         isotope_evidence_score            eips_evidence_score 
#>                           0.20                           0.15 
#> standard_adduct_recovery_score         family_coherence_score 
#>                           0.15                           0.10 
```
