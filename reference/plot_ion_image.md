# Plot an ion image

Plots the spatial ion image (per-pixel intensity) of a single feature in
a peak matrix, using `ggplot2`.

## Usage

``` r
plot_ion_image(
  pkm,
  idx = NULL,
  mz = NULL,
  tol_ppm = 10,
  palette = "viridis",
  title = NULL,
  flip_y = TRUE,
  clip_quantile = NULL,
  show_legend = TRUE,
  group = NULL
)
```

## Arguments

- pkm:

  Peak matrix list containing at least `mass` (numeric vector),
  `intensity` (numeric pixel x feature matrix) and `pos` (matrix or
  data.frame with `x`, `y` pixel coordinate columns).

- idx:

  Optional integer. Column index of the feature to plot in
  `pkm$mass`/`pkm$intensity` (the same indexing convention used by
  `feature_summary$idx`, `morph_results$idx_M0`/`idx_cand`,
  `adduct_edges$idx_i`/`idx_j`, `candidate_annotations$feature_idx` and
  `adduct_fam$family_members$idx`). This is the preferred way to select
  a feature. Exactly one of `idx` or `mz` must be supplied.

- mz:

  Optional numeric. Target m/z used to look up the closest feature in
  `pkm$mass` within `tol_ppm`, for exploratory use when a feature index
  is not already known. Exactly one of `idx` or `mz` must be supplied.

- tol_ppm:

  PPM tolerance used when searching by `mz`. Ignored (without a warning)
  when `idx` is supplied.

- palette:

  Colour palette for the intensity scale. `"viridis"` (default) uses the
  blue-to-yellow
  [`ggplot2::scale_fill_viridis_c()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
  scale. Any other value is passed through as the `option` argument of
  `scale_fill_viridis_c()` (e.g. `"magma"`, `"inferno"`).

- title:

  Optional plot title. If `NULL` (default), a title is generated from
  the plotted feature's m/z - when the feature was located via `mz`, the
  title reflects the m/z actually found, not the one requested.

- flip_y:

  Logical. If `TRUE` (default), the pixel y coordinate is negated before
  plotting, matching the top-left image origin convention used elsewhere
  in the package.

- clip_quantile:

  Optional numeric in `(0, 1]`. If supplied, intensity values above this
  quantile are clipped before mapping to colour, so a handful of
  saturated pixels do not flatten the rest of the colour scale.

- show_legend:

  Logical. If `FALSE`, the intensity colour legend is omitted.

- group:

  Optional grouping of the pixels into separate images, used to draw one
  panel per image side by side with
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  If `NULL` (default), the image identity stored in `pkm$image_id`
  (added by
  [`cardinal_to_peakmatrix()`](https://nadiamolto.github.io/PeakGuideR/reference/cardinal_to_peakmatrix.md))
  or described by `pkm$numPixels`/`pkm$names` (native rMSI2 peak
  matrices) is used when it contains more than one image; otherwise a
  single panel is drawn. Alternatively, a vector with one value per
  pixel (`nrow(pkm$intensity)`), or `FALSE` to always draw a single
  panel. Pixel coordinates are not required to be unique across images,
  so overlapping coordinates between images are drawn in separate panels
  instead of on top of each other.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_pkm, package = "PeakGuideR")

plot_ion_image(example_pkm, idx = 1)
plot_ion_image(example_pkm, mz = 762.6, tol_ppm = 50)
} # }
```
