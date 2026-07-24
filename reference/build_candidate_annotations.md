# Build the unified long-format candidate annotation table

Builds one long candidate-annotation table combining adduct-family
evidence and single-adduct hypotheses (see
[`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)),
with database and standard-adduct identity resolved via
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)
instead of the identifier-hierarchy-free join used previously. This is
the table returned as `candidate_annotations` by
[`run_peakguider_workflow()`](https://nadiamolto.github.io/PeakGuideR/reference/run_peakguider_workflow.md).
[`build_neutral_mass_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_neutral_mass_candidates.md)
is a per-`neutral_mass_id` summary view derived from this table.

Each row represents one neutral mass, candidate identity (from
`compound_db`, `standards_db`, or both fused via a strong identifier
match) and evidence source. A match on molecular formula alone is never
fused into a single identity; see
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md).

## Usage

``` r
build_candidate_annotations(
  adduct_fam,
  feature_summary,
  pkm = NULL,
  ion_mode = c("pos", "neg"),
  matrix = NULL,
  adducts = NULL,
  compound_db = NULL,
  standards_db = NULL,
  ppm_tol = 5,
  neutral_cluster_ppm = 5,
  top_n = 10L,
  weights = default_priority_weights(),
  recovery_threshold = 0.5,
  adduct_min_score = 0.5,
  ambiguity_gap = 0.05,
  include_single_adduct = TRUE,
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

- pkm:

  Optional peak matrix (`mass`, `intensity`). Required to include
  single-adduct hypotheses (see `include_single_adduct`); if `NULL`, the
  table contains family-derived candidates only.

- ion_mode:

  `"pos"` or `"neg"`.

- matrix:

  Matrix name. Standard-adduct support is applied only when
  `matrix = "HCCA"`. Use `NULL` to skip it.

- adducts:

  Optional adduct definition data.frame, used only for single-adduct
  hypotheses. If `NULL`, `default_adducts(ion_mode)` is used.

- compound_db:

  Compound mass database. If `NULL`, the included example database is
  loaded.

- standards_db:

  Standard adduct library. If `NULL` and `matrix = "HCCA"`, the included
  example library is loaded.

- ppm_tol:

  PPM tolerance for compound/standard mass matching.

- neutral_cluster_ppm:

  PPM tolerance used to cluster inferred neutral masses.

- top_n:

  Maximum number of compound candidates kept per neutral mass. Use
  `NULL` to retain all candidates within `ppm_tol`.

- weights:

  Named numeric vector of weights for `priority_score`; see
  [`default_priority_weights()`](https://nadiamolto.github.io/PeakGuideR/reference/default_priority_weights.md).

- recovery_threshold:

  Minimum `standard_adduct_recovery_score` (`NA` treated as `0`)
  required for
  `confidence_class = "identity_confirmed_adduct_recovered"` when
  `source == "both"`.

- adduct_min_score:

  Minimum `adduct_spatial_score` (`NA` treated as `0`) counted as
  multi-evidence support for `broad_db_only` candidates. Reuses the same
  threshold as the `adduct_min_score` workflow parameter.

- ambiguity_gap:

  Minimum `priority_score` gap between the top-1 and top-2 candidate of
  the same `neutral_mass_id` required to *not* flag
  `ambiguous_isomeric`.

- include_single_adduct:

  Logical. If `FALSE`, or if `pkm` is `NULL`, single-adduct hypotheses
  are not computed.

- quiet:

  Logical. If `FALSE`, database loading functions may print notices.

## Value

A data.frame with one row per neutral mass, feature(s) and candidate
identity, including `neutral_mass_id`, `neutral_mass_consensus`,
`feature_idx`, `inferred_adduct`, `source`, `identity_match_type`,
`possible_link_formula_only`, `hypothesis_origin`, `broad_db_*`,
`standard_db_*`, `mass_error_score`, `adduct_spatial_score`,
`isotope_evidence_score`, `eips_evidence_score`,
`standard_adduct_recovery_score`, `family_coherence_score`,
`priority_score`, `confidence_class` and `ambiguous_isomeric`.

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

candidate_annotations <- build_candidate_annotations(
  adduct_fam = adduct_fam,
  feature_summary = feature_summary,
  pkm = example_pkm,
  ion_mode = "pos",
  matrix = "HCCA"
)
head(candidate_annotations)
} # }
```
