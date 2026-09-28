# Derive a per-pixel image identity from a peak matrix

Internal helper returning a factor with one entry per pixel (row of
`pkm$intensity`) that identifies which image/run each pixel belongs to,
so that peak matrices coming from different input paths share one
representation.

The identity is taken, in order of preference, from:

1.  `pkm$image_id`, as added by
    [`cardinal_to_peakmatrix()`](https://nadiamolto.github.io/PeakGuideR/reference/cardinal_to_peakmatrix.md).

2.  `pkm$numPixels` and `pkm$names`, as carried natively by rMSI2
    `rMSIprocPeakMatrix` objects with more than one image: pixels are
    concatenated in contiguous row-blocks following that order.

If neither describes more than one image, every pixel is assigned to a
single image (named after `pkm$names[1]` when available).

## Usage

``` r
derive_image_id_from_pkm(pkm)
```

## Arguments

- pkm:

  Peak matrix list containing at least `intensity`.

## Value

A factor of length `nrow(pkm$intensity)`, with levels in order of first
appearance.
