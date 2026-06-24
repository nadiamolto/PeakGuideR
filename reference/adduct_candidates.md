# Detects adduct pairwise candidates

Detects pairwise adduct candidates by combining:

- expected m/z differences between adduct species

- spatial/intensity similarity between features

- agreement between inferred neutral masses.

The function returns candidate adduct edges that can later be grouped
into larger ion families or pseudo-compounds.

## Usage

``` r
adduct_candidates(
  pkm,
  ion_mode = c("pos", "neg"),
  adducts = NULL,
  feature_idx = NULL,
  tol_ppm = 5,
  neutral_tol_ppm = 5,
  method = c("pearson", "cosine", "spearman"),
  transform = c("none", "log1p", "zscore"),
  min_quantile = 0.01,
  clip_negatives = TRUE,
  min_score_spatial = 0.5
)
```

## Arguments

- pkm:

  A list with:

  - mass: numeric vector of m/z values.

  - intensity: numeric matrix pixels by features.

- ion_mode:

  Ionization mode. Either "pos" or "neg".

- adducts:

  Optional data-frame of adduct definitions. If NULL, default_adducts()
  is used (adducts usually found in MALDI-MSI).

- feature_idx:

  Optional integer vector of feature indices to evaluate. If NULL, all
  features are used.

- tol_ppm:

  Numeric. PPM tolerance used to search candidate peak pairs. Default 5.

- neutral_tol_ppm:

  Numeric. PPM tolerance used to score agreement between inferred
  neutral masses. Default `5`. This value was chosen for high-resolution
  Orbitrap data and should be adjusted for other instruments or
  preprocessing settings.

- method:

  Similarity metric between feature intensity profiles: "pearson",
  "cosine" or "spearman". Default "pearson".

- transform:

  Intensity transformation before similarity calculation: "none",
  "log1p" or "zscore". Default "none".

- min_quantile:

  Numeric in the range 0 to 1. Feature-wise low-intensity quantile
  filter applied after retaining co-detected pixels. Pixels are kept
  only when both features are above their respective threshold. Default
  `0.01`.

- clip_negatives:

  Logical. If TRUE, negative intensities are truncated to zero before
  transformation. Default TRUE.

- min_score_spatial:

  Numeric in the range 0 to 1. Minimum spatial score required to keep an
  adduct edge. Default 0.5.

## Value

A data.frame with one row per candidate adduct relation:

- idx_i, mz_i: first feature index and m/z

- idx_j, mz_j: second feature index and m/z

- adduct_i, adduct_j: adduct hypotheses

- delta_theo, delta_obs, delta_err_ppm

- neutral_mass_i, neutral_mass_j, neutral_mass_mean

- neutral_err_ppm

- score_spatial, score_mass, score_adduct

- is_valid_adduct

## Details

The function assumes that the mass column in the adduct table stores the
net m/z shift relative to the neutral mass:

\$\$mz = M + shift\$\$

so the neutral mass is recovered as:

\$\$M = mz - shift\$\$

For each pair of adducts within the selected ionization mode, the
expected m/z difference is:

\$\$\Delta m/z = shift_b - shift_a\$\$

Candidate feature pairs matching that difference are then evaluated
using spatial/intensity similarity and neutral-mass consistency.

The reported `score_adduct` corresponds to the spatial similarity score
after filtering candidate pairs by the expected adduct m/z difference.
The neutral-mass consistency score is returned as `score_mass` for
inspection and downstream filtering.

## Examples

``` r
if (FALSE) { # \dontrun{
adduct_res <- adduct_candidates(
  pkm,
  ion_mode = "pos"
)

extra_adducts <- rbind(
  default_adducts(),
  data.frame(
    name = "[M-H]-",
    mode = "neg",
    mass = -1.007276,
    stringsAsFactors = FALSE
  )
)

adduct_res_neg <- adduct_candidates(
  pkm,
  ion_mode = "neg",
  adducts = extra_adducts
)
} # }
```
