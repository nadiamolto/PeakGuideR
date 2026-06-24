# Detect isotope candidates using mass differences and ion-image morphology

Detects potential isotopic relationships between peaks in a Mass
Spectrometry Imaging (MSI) peak matrix based on expected mass
differences and spatial colocalization (correlation). It computes
correlation-like morphology metrics (Pearson, cosine, or Spearman),
optionally evaluates spatial consistency across tiles, and returns all
isotopic pairs that satisfy mass tolerance and a minimum similarity
threshold.

## Usage

``` r
iso_morphology_candidates(
  pkm,
  prefer_mode = c("ppm", "dp"),
  tol_ppm = 5,
  tol_dp = 3L,
  method = c("pearson", "cosine", "spearman"),
  transform = c("none", "log1p", "zscore"),
  min_quantile = 0.01,
  clip_negatives = TRUE,
  use_tiles = TRUE,
  tiles = 9L,
  tile_min_pixels = 50L,
  tile_blend = c("median", "p25", "pass_rate"),
  tile_alpha = 0.8,
  tile_threshold = 0.8,
  min_score_keep = 0.2
)
```

## Arguments

- pkm:

  A list with:

  - mass: numeric vector of m/z values.

  - intensity: matrix pixel x feature.

  - optional pos: data.frame/matrix with x and y coordinates for pixels
    (without it, tile strategy is disabled with an info message).
    Position information is recommended to get more reliable results.

- prefer_mode:

  "ppm" (default) or "dp". Tolerance mode.

- tol_ppm:

  Numeric, mass window (±ppm). Default 5.

- tol_dp:

  Integer, window in indices if prefer_mode="dp". Default 3L.

- method:

  Character, morphology metric: "pearson" (default), "cosine", or
  "spearman".

- transform:

  One of "none", "log1p", "zscore". Default "none".

- min_quantile:

  Numeric in the range 0 to 1, pooled lower quantile filter. Default
  0.01.

- clip_negatives:

  Logical, clip negative values before transform. Default TRUE.

- use_tiles:

  Logical, enable spatial tile evaluation if pkm\$pos exists. Default
  TRUE.

- tiles:

  Integer, approximate number of subregions. Default 9L.

- tile_min_pixels:

  Integer, minimum pixels per tile. Default 50L.

- tile_blend:

  One of "median", "p25", or "pass_rate". Default "median".

- tile_alpha:

  Numeric in the range 0 to 1, blending between global and tile score.
  Default 0.8.

- tile_threshold:

  Numeric in the range 0 to 1, threshold for "pass_rate". Default 0.8.

- min_score_keep:

  Numeric in the range 0 to 1, score threshold. Default 0.2.

## Value

A data.frame with:

- idx_M0, mz_M0: anchor index and m/z.

- iso_type: isotope type (C13_M1, C13_M2, S34, Cl37, Br81, O18,15N).

- element: (C.N,O,S,Cl,Br)

- k: isotope level (1 or 2).

- z: charge state (always 1).

- idx_cand, mz_cand: candidate index and m/z.

- score_global, tile_summary, tile_sd, tile_consistency, score_final.

- mass_err_da, mass_err_ppm: candidate error vs theoretical mass.

- mass dev_err, mass_dev_score: absolute mass error and corresponding
  mass deviation score.

## Details

Charge state is **fixed to z = 1** throughout the method (no user
control).

For each monoisotopic feature (M+0), the function searches isotopic
candidates at theoretical mass differences for z = 1:

- **13C series:** M+1 = 1.003355, M+2 = 2 × 1.003355

- **Other isotopes (M+X level):** 34S = 1.99580, 37Cl = 1.99705, 81Br =
  1.99795, 18O = 2.004245

**Tolerance modes:**

- "ppm" (default): find candidates within ±tol_ppm ppm around the
  theoretical mass.

- "dp": find candidates within ±tol_dp data points in the ordered mass
  axis.

**Preprocessing steps per intensity pair:**

1.  Keep pixels where both features are detected (x!=0 & y!=0) \#' 2.
    Apply feature-wise low-intensity quantile filters to the retained
    pixels. Pixels are kept only when both features are above their
    respective threshold.

2.  Optionally clip negatives (clip_negatives=TRUE).

3.  Apply intensity transform (none, log1p, or zscore). Scoring metrics
    should be considered before setting the transformation method:
    zscore is recommended for pearson correlation, log1p or none for
    cosine similarity and none for spearman. Normalization can be
    applied (e.g. TIC) before it but should be considered regarding
    result interpretation.

**Scoring metrics (method):**

- "pearson": squared Pearson correlation (R²), 0 to 1.

- "cosine": cosine similarity.

- "spearman": squared Spearman rank correlation (negatives set to 0).

**Tile mode:** if use_tiles=TRUE and pixel coordinates are present
(pkm\$pos\$x,y), the sample is split into approximately
sqrt(tiles)×sqrt(tiles) subregions. For each tile, the same correlation
is computed, summarizing via:

- "median": median of tile scores,

- "p25": 25th percentile (robust lower-bound),

- "pass_rate": fraction of tiles with score ≥ tile_threshold.

A consistency factor is computed as: \$\$ tile\\consistency = \max(0,
1 - sd(tile\\scores)) \$\$

The final score blends the global morphology score with the
tile-consistency factor:

\$\$score\\final = \alpha \cdot score\\global + (1-\alpha) \cdot
tile\\consistency\$\$

The `tile_alpha` parameter controls the contribution of the global score
relative to the tile-consistency component. The selected `tile_blend`
summary is returned as `tile_summary` for inspection.

**Output policy:** Returns *all* isotopic candidate pairs (M+0... M+k)
that fall within the mass tolerance and have score_final ≥
min_score_keep. A dataframe with as many rows as isotopes are found.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- iso_morphology_candidates(
  pkm,
  prefer_mode   = "ppm",
  tol_ppm       = 5,
  method        = "pearson",
  transform     = "log1p",
  use_tiles     = TRUE,
  tiles         = 9L,
  tile_blend    = "median",
  tile_alpha    = 0.8,
  min_score_keep= 0.2
)
head(res)
} # }
```
