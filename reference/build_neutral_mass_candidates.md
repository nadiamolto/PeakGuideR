# Build neutral-mass candidate table

Builds a neutral-mass candidate table from adduct families.

This function is a per-`neutral_mass_id` summary view derived from
[`build_candidate_annotations()`](https://nadiamolto.github.io/PeakGuideR/reference/build_candidate_annotations.md)
(family-derived candidates only; see `include_single_adduct` there for
the full single-adduct-aware table returned as `candidate_annotations`
by
[`run_peakguider_workflow()`](https://nadiamolto.github.io/PeakGuideR/reference/run_peakguider_workflow.md)).
Compound and standard-adduct identity are resolved with
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md),
so a shared molecular formula alone never fuses two different compounds
into the same `has_standard_compound_match` (see the `dontrun` example
in
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)
for details).

The output contains one row per inferred neutral mass and
compound/standard candidate. Candidate compounds are putative mass
matches, not definitive identifications.

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
  candidate_annotations = NULL,
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

  Maximum number of compound candidates kept per neutral mass. Use
  `NULL` to retain all candidates within `ppm_tol`.

- candidate_annotations:

  Optional pre-computed output from
  [`build_candidate_annotations()`](https://nadiamolto.github.io/PeakGuideR/reference/build_candidate_annotations.md).
  If `NULL` (the default), the function computes its own family-only
  candidate table internally, exactly as before. If supplied (for
  example the same table already computed by
  [`run_peakguider_workflow()`](https://nadiamolto.github.io/PeakGuideR/reference/run_peakguider_workflow.md),
  possibly including single-adduct rows), it is reused instead of
  recomputing the compound/standard matching and identity resolution
  from scratch; rows other than `hypothesis_origin == "family"` are
  filtered out first, so the returned table is unaffected by whether
  single-adduct hypotheses were included upstream.

- quiet:

  Logical. If `FALSE`, database loading functions may print notices.

## Value

A data.frame with one row per neutral mass and compound/standard
candidate.

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

neutral_mass_candidates <- build_neutral_mass_candidates(
  adduct_fam = adduct_fam,
  feature_summary = feature_summary,
  ion_mode = "pos",
  matrix = "HCCA"
)
head(neutral_mass_candidates)
} # }
```
