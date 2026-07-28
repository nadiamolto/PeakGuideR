# Plot an isotope pair (and, when available, its M+2 satellite)

Plots the ion image of a candidate monoisotopic feature (M0) side by
side with its candidate isotope satellite (M+k), reusing
[`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md)
for each panel. If a
[`cir_score()`](https://nadiamolto.github.io/PeakGuideR/reference/cir_score.md)
row is supplied and it has a supported M+2 satellite
(`has_C13_M2 == TRUE`), a third panel for that M+2 feature is added
automatically.

## Usage

``` r
plot_isotope_pair(
  pkm,
  idx_M0,
  idx_cand,
  morph_row = NULL,
  cir_row = NULL,
  eips_row = NULL,
  ...
)
```

## Arguments

- pkm:

  Peak matrix list, see
  [`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md).

- idx_M0:

  Integer. Column index of the monoisotopic (M0) feature.

- idx_cand:

  Integer. Column index of the candidate isotope satellite feature
  (M+k).

- morph_row:

  Optional single-row data.frame - one row of the output of
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md) -
  used (if `cir_row` and `eips_row` are not supplied) to annotate the
  plot subtitle with `score_final`.

- cir_row:

  Optional single-row data.frame - one row of the output of
  [`cir_score()`](https://nadiamolto.github.io/PeakGuideR/reference/cir_score.md).
  When supplied, its `R_obs`/`R_theo`/`cir_class` are used for the plot
  subtitle (taking priority over `eips_row`/`morph_row`), and an M+2
  panel is added automatically when `has_C13_M2` is `TRUE`.

- eips_row:

  Optional single-row data.frame - one row of the output of
  [`eips_score()`](https://nadiamolto.github.io/PeakGuideR/reference/eips_score.md).
  Used for the plot subtitle when `cir_row` is not supplied.

- ...:

  Additional arguments passed through to
  [`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md)
  (e.g. `palette`, `clip_quantile`, `flip_y`, `show_legend`).

## Value

A `patchwork` object combining the M0/M+k (and, when available, M+2) ion
image panels.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_pkm, package = "PeakGuideR")

morph_results <- iso_morphology_candidates(example_pkm, prefer_mode = "ppm")
cir_results <- cir_score(morph_results, example_pkm)

one_row <- cir_results[1, ]
plot_isotope_pair(
  example_pkm,
  idx_M0 = one_row$idx_M0,
  idx_cand = one_row$idx_M1,
  cir_row = one_row
)
} # }
```
