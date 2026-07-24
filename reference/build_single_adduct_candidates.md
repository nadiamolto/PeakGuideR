# Build single-adduct candidate hypotheses for unassigned features

Builds candidate neutral-mass hypotheses for m/z features that were not
assigned to any adduct family by
[`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md),
excluding features already explained as isotope satellite peaks
(`is_c13_m1`, `is_c13_m2`, `is_eips_isotope` in `feature_summary`).

For each unassigned feature and each adduct compatible with `ion_mode`,
a neutral-mass hypothesis is inferred
(`neutral_mass = mz - adduct_mass`, the same convention as elsewhere in
PeakGuideR). These hypotheses are then clustered by ppm tolerance using
the same mechanism as
[`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md)
/
[`build_neutral_mass_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_neutral_mass_candidates.md)
(see
[`assign_neutral_mass_clusters()`](https://nadiamolto.github.io/PeakGuideR/reference/assign_neutral_mass_clusters.md)),
and share the **same** `neutral_mass_id` space as adduct families: a
hypothesis that falls within `neutral_cluster_ppm` of an
already-detected family's neutral mass is folded into that family's
`neutral_mass_id` (`hypothesis_origin = "single_adduct_clustered"`)
instead of becoming a separate entry. Isolated hypotheses that do not
fall into any existing family are clustered only among themselves and
receive new `neutral_mass_id` values
(`hypothesis_origin = "single_adduct_isolated"`).

An `adduct_spatial_score` is computed for every hypothesis by
correlating the feature's intensity image against every other feature
already linked to the same `neutral_mass_id` (family members and/or
other single-adduct hypotheses), reusing the spatial-similarity
utilities in `R/spatial_similarity.R`. A hypothesis with no such partner
(a truly isolated cluster of size one) gets `adduct_spatial_score = NA`,
not `0`.

Candidates are then retrieved by mass from `compound_db`
([`match_compound_db_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_compound_db_by_mass.md))
and, when `matrix = "HCCA"`, from `standards_db`
([`match_standards_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_standards_by_mass.md)),
and identity is resolved with
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)
exactly as for family-derived candidates.
`standard_adduct_recovery_score` is computed against the
detected-adducts set enriched by the neutral-mass clustering above
(family members' adducts and every single-adduct hypothesis sharing the
same `neutral_mass_id`).

## Usage

``` r
build_single_adduct_candidates(
  pkm,
  feature_summary,
  adduct_fam,
  ion_mode = c("pos", "neg"),
  matrix = NULL,
  adducts = NULL,
  compound_db = NULL,
  standards_db = NULL,
  ppm_tol = 5,
  neutral_cluster_ppm = 5,
  top_n = 10L,
  method = c("pearson", "cosine", "spearman"),
  transform = c("none", "log1p", "zscore"),
  min_quantile = 0.01,
  clip_negatives = TRUE,
  quiet = FALSE
)
```

## Arguments

- pkm:

  A peak matrix list with `mass` and `intensity`.

- feature_summary:

  Output from
  [`build_feature_summary()`](https://nadiamolto.github.io/PeakGuideR/reference/build_feature_summary.md).

- adduct_fam:

  Output from
  [`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md).

- ion_mode:

  `"pos"` or `"neg"`.

- matrix:

  Matrix name. Standard-adduct support is applied only when
  `matrix = "HCCA"`. Use `NULL` to skip it.

- adducts:

  Optional adduct definition data.frame. If `NULL`,
  `default_adducts(ion_mode)` is used.

- compound_db:

  Compound mass database. If `NULL`, the included example database is
  loaded.

- standards_db:

  Standard adduct library. If `NULL` and `matrix = "HCCA"`, the included
  example library is loaded.

- ppm_tol:

  PPM tolerance for compound/standard mass matching.

- neutral_cluster_ppm:

  PPM tolerance used to cluster inferred neutral masses, including
  folding hypotheses into existing adduct families.

- top_n:

  Maximum number of compound candidates kept per neutral mass. Use
  `NULL` to retain all candidates within `ppm_tol`.

- method:

  Spatial similarity metric: `"pearson"`, `"cosine"` or `"spearman"`.

- transform:

  Intensity transformation before similarity calculation: `"none"`,
  `"log1p"` or `"zscore"`.

- min_quantile:

  Numeric in the range 0 to 1. Feature-wise low-intensity quantile
  filter used by the spatial similarity score.

- clip_negatives:

  Logical. If `TRUE`, negative intensities are truncated to zero before
  transformation.

- quiet:

  Logical. If `FALSE`, database loading functions may print notices.

## Value

A data.frame with one row per feature, inferred adduct and candidate
identity: `neutral_mass_id`, `neutral_mass_consensus`, `feature_idx`,
`inferred_adduct`, `hypothesis_origin`, `adduct_spatial_score`,
`source`, `identity_match_type`, `possible_link_formula_only`,
`candidate_ppm_error`, `candidate_neutral_mass`, `broad_db_*`,
`standard_db_*` and `standard_adduct_recovery_score`.

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
feature_summary <- build_feature_summary(
  relation_table = relation_table,
  adduct_fam = adduct_fam,
  pkm = example_pkm
)

single_adduct_hyp <- build_single_adduct_candidates(
  pkm = example_pkm,
  feature_summary = feature_summary,
  adduct_fam = adduct_fam,
  ion_mode = "pos",
  matrix = "HCCA"
)
head(single_adduct_hyp)
} # }
```
