# Build a unified peak relation table

Builds a unified table of peak-to-peak relationships detected by
PeakGuideR.

This table integrates internal evidence from C13 isotope validation, C13
M+2 support, elemental isotope-pattern support and adduct-family
relationships.

Each relation keeps the original score from the evidence layer that
generated it. Therefore, `evidence_score` should be interpreted together
with `score_type`.

## Usage

``` r
build_relation_table(
  cir_results = NULL,
  eips_results = NULL,
  adduct_fam = NULL,
  only_valid = TRUE
)
```

## Arguments

- cir_results:

  Output from
  [`cir_score()`](https://nadiamolto.github.io/PeakGuideR/reference/cir_score.md).

- eips_results:

  Output from
  [`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md).

- adduct_fam:

  Output from
  [`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md).

- only_valid:

  Logical. If `TRUE`, keeps only valid relations.

## Value

A data.frame with one row per peak-to-peak relation.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_pkm, package = "PeakGuideR")

morph_results <- iso_morphology_candidates(example_pkm, prefer_mode = "ppm")
cir_results <- cir_score(morph_results, example_pkm)
eips_results <- eips_score(
  morph_results, example_pkm,
  ion_mode = "pos", cir_df = cir_results, morph_df = morph_results
)

adduct_edges <- adduct_candidates(example_pkm, ion_mode = "pos")
adduct_fam <- adduct_families(adduct_edges)

relation_table <- build_relation_table(
  cir_results = cir_results,
  eips_results = eips_results,
  adduct_fam = adduct_fam
)
head(relation_table)
} # }
```
