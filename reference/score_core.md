# Spatial/intensity similarity score in the range 0 to 1

Shared internal similarity metric used by
[`adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_candidates.md)
and
[`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)
on preprocessed intensity vectors (see
[`preprocess_xy()`](https://nadiamolto.github.io/PeakGuideR/reference/preprocess_xy.md)).

## Usage

``` r
score_core(x, y, method = c("pearson", "cosine", "spearman"))
```

## Arguments

- x, y:

  Numeric vectors, typically the output of
  [`preprocess_xy()`](https://nadiamolto.github.io/PeakGuideR/reference/preprocess_xy.md).

- method:

  Similarity metric: `"pearson"`, `"cosine"` or `"spearman"`.

## Value

A single numeric value in the range 0 to 1, or `NA_real_`.
