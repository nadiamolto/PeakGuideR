# Run the PeakGuideR annotation workflow

Runs the main PeakGuideR workflow from either a PeakGuideR peak matrix
object or a Cardinal MSI object.

The workflow includes isotope morphology detection, carbon isotope-ratio
validation, elemental isotope-pattern support, adduct candidate
detection, adduct-family grouping, relation-table construction,
feature-level summarization and neutral-mass candidate matching.

## Usage

``` r
run_peakguider_workflow(
  pkm,
  ion_mode = c("pos", "neg"),
  matrix = NULL,
  adducts = NULL,
  compound_db = NULL,
  standards_db = NULL,
  morph_prefer_mode = c("ppm", "dp"),
  morph_method = c("pearson", "cosine", "spearman"),
  morph_transform = c("none", "log1p", "zscore"),
  morph_tile_blend = c("median", "p25", "pass_rate"),
  iso_min_score = 0.6,
  cir_rel_tol = 0.3,
  eips_rel_tol = 0.3,
  ratio_method = c("sum", "mean", "median"),
  adduct_tol_ppm = 5,
  adduct_neutral_tol_ppm = 5,
  adduct_method = c("pearson", "cosine", "spearman"),
  adduct_transform = c("none", "log1p", "zscore"),
  adduct_min_quantile = 0.01,
  adduct_clip_negatives = TRUE,
  adduct_min_score = 0.5,
  neutral_cluster_ppm = 5,
  candidate_ppm_tol = 5,
  top_n = 10L,
  include_single_adduct = TRUE,
  quiet = FALSE,
  multi_image = FALSE
)
```

## Arguments

- pkm:

  rMSI2 peak matrix object, or a supported Cardinal
  `MSImagingExperiment` object. Cardinal objects are converted
  internally using
  [`cardinal_to_peakmatrix()`](https://nadiamolto.github.io/PeakGuideR/reference/cardinal_to_peakmatrix.md).

- ion_mode:

  Ion mode, either `"pos"` or `"neg"`.

- matrix:

  Matrix name. Use `"HCCA"` to enable HCCA-specific standard adduct
  support in neutral-mass candidate matching.

- adducts:

  Optional adduct definition table. If `NULL`, PeakGuideR uses
  `default_adducts(ion_mode)`. Users can inspect and modify the default
  adduct table with
  [`default_adducts()`](https://nadiamolto.github.io/PeakGuideR/reference/default_adducts.md).

- compound_db:

  Optional compound mass database.

- standards_db:

  Optional standard adduct library.

- morph_prefer_mode:

  Mass-deviation preference mode used by
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md),
  either `"ppm"` or `"dp"`.

- morph_method:

  Spatial similarity method used by
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md).

- morph_transform:

  Intensity transformation used by
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md).

- morph_tile_blend:

  Method used to combine tile-level morphology scores.

- iso_min_score:

  Minimum isotope morphology score used for CIR/EIPS input.

- cir_rel_tol:

  Relative tolerance for carbon isotope-ratio validation.

- eips_rel_tol:

  Relative tolerance for EIPS validation.

- ratio_method:

  Ratio aggregation method.

- adduct_tol_ppm:

  PPM tolerance for adduct candidate detection.

- adduct_neutral_tol_ppm:

  PPM tolerance for adduct neutral-mass consistency.

- adduct_method:

  Spatial similarity method used for adduct detection.

- adduct_transform:

  Intensity transformation used for adduct detection.

- adduct_min_quantile:

  Minimum quantile used for adduct spatial vectors.

- adduct_clip_negatives:

  Logical. If `TRUE`, negative transformed values are clipped in adduct
  detection.

- adduct_min_score:

  Minimum spatial score for adduct candidates/families.

- neutral_cluster_ppm:

  PPM tolerance used to cluster inferred neutral masses.

- candidate_ppm_tol:

  PPM tolerance for compound candidate matching.

- top_n:

  Maximum number of compound candidates per neutral mass. Use `NULL` to
  retain all candidates within `candidate_ppm_tol`.

- include_single_adduct:

  Logical. If `TRUE` (the default), `candidate_annotations` also
  includes single-adduct hypotheses for features not assigned to any
  adduct family (see
  [`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)).
  If `FALSE`, only family-derived candidates are included.

- quiet:

  Logical. If `FALSE`, prints progress messages.

- multi_image:

  Logical. Set to `TRUE` when `pkm` contains more than one image/run
  concatenated by pixel (for example a Cardinal object with several
  runs, or an rMSI2 peak matrix with several images). Annotation stays
  fully joint: correlation, EIPS, adduct detection and candidate scoring
  pool all pixels of all images, and a single `candidate_annotations`
  table is returned. Because pixel coordinates are local to each image
  and may overlap between images, **the tile-based spatial-consistency
  step of
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)
  is disabled** (`use_tiles = FALSE`), so
  `morph_results$tile_consistency` is `NA` and the isotope-morphology
  score relies on the global spatial score only. The per-pixel image
  identity is attached to the returned `pkm` as `image_id` so that
  [`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md)
  shows each image in its own panel. The default `FALSE` leaves the
  workflow unchanged.

## Value

A list with all main PeakGuideR workflow outputs:

- `morph_results`: isotope-morphology candidates (the output of
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)),
  one candidate peak pair per row.

- `cir_results`: carbon isotope-ratio (CIR) validation for the M0/M+1
  pairs in `morph_results`, including optional M+2 support.

- `eips_results`: elemental isotope-pattern (EIPS) validation for N, O,
  S, Cl and Br.

- `adduct_edges`: candidate adduct peak pairs, before grouping into
  families.

- `adduct_families`: adduct pairs grouped into families by consensus
  neutral mass.

- `relation_table`: a single table unifying every peak-to-peak relation
  (CIR, isotope morphology, EIPS, adduct family).

- `feature_summary`: one row per detected peak/feature, with the role it
  plays in each type of evidence.

- `neutral_mass_candidates`: a per-inferred-neutral-mass summary.

- `candidate_annotations`: the final compound candidate annotations (one
  row per neutral mass and candidate identity).

- `pkm`: the peak matrix object used for the analysis (with `image_id`
  added when `multi_image = TRUE`) - see the "Memory usage" section
  below for why it can be worth freeing.

- `parameters`: a list with every parameter value used for the call.

## Memory usage

The returned `pkm` element is a copy of the `pkm` argument (with
`image_id` attached when `multi_image = TRUE`), not a reference to it.
For a real MSI dataset of several gigabytes, keeping both the original
`pkm` object and `res$pkm` alive at the same time doubles that memory
footprint for no benefit. Once `res$pkm` is no longer needed from `res`
itself, it can be dropped with `res$pkm <- NULL` (followed by
[`gc()`](https://rdrr.io/r/base/gc.html) if the memory needs to be
returned to the OS immediately rather than whenever R next collects);
the original `pkm` object already in the user's environment can still be
passed directly to
[`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md),
[`plot_isotope_pair()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_isotope_pair.md)
and
[`plot_adduct_family()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_adduct_family.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_pkm, package = "PeakGuideR")

res <- run_peakguider_workflow(
  pkm = example_pkm,
  ion_mode = "pos",
  matrix = "HCCA"
)

names(res)
head(res$candidate_annotations)
} # }
```
