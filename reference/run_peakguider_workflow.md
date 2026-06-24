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
  quiet = FALSE
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

  Maximum number of compound candidates per neutral mass.

- quiet:

  Logical. If `FALSE`, prints progress messages.

## Value

A list with all main PeakGuideR workflow outputs.
