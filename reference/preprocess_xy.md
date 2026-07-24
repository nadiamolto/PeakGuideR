# Preprocess two intensity vectors for spatial similarity scoring

Shared internal preprocessing used by
[`adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_candidates.md)
and
[`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)
to prepare a pair of feature intensity vectors before computing a
spatial/intensity similarity score.

## Usage

``` r
preprocess_xy(
  x,
  y,
  min_quantile = 0.01,
  clip_negatives = TRUE,
  transform = c("none", "log1p", "zscore")
)
```

## Arguments

- x, y:

  Numeric intensity vectors (same length).

- min_quantile:

  Numeric in the range 0 to 1. Feature-wise low-intensity quantile
  filter applied after retaining co-detected pixels.

- clip_negatives:

  Logical. If `TRUE`, negative intensities are truncated to zero before
  transformation.

- transform:

  Intensity transformation: `"none"`, `"log1p"` or `"zscore"`.

## Value

A list with `x` and `y` (preprocessed vectors), or `NULL` if fewer than
3 valid pixels remain.
