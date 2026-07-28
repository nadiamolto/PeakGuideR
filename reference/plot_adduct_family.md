# Plot all members of an adduct family

Plots the ion image of every feature belonging to a given adduct family
(as returned by
[`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md)),
reusing
[`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md)
for each panel and combining them with `patchwork`.

## Usage

``` r
plot_adduct_family(pkm, family_id, adduct_fam, ncol = NULL, ...)
```

## Arguments

- pkm:

  Peak matrix list, see
  [`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md).

- family_id:

  Adduct family ID to plot (matches
  `adduct_fam$family_members$family_id` /
  `adduct_fam$family_summary$family_id`).

- adduct_fam:

  List returned by
  [`adduct_families()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_families.md),
  containing `family_members` and `family_summary`.

- ncol:

  Optional integer. Number of columns used to arrange the panels. If
  `NULL` (default), a reasonable number of columns is chosen
  automatically from the family size (`ceiling(sqrt(n))`).

- ...:

  Additional arguments passed through to
  [`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md)
  (e.g. `palette`, `clip_quantile`, `flip_y`, `show_legend`).

## Value

A `patchwork` object combining one ion image panel per family member.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_pkm, package = "PeakGuideR")

adduct_res <- adduct_candidates(example_pkm, ion_mode = "pos")
adduct_fam <- adduct_families(adduct_res)

plot_adduct_family(
  example_pkm,
  family_id = adduct_fam$family_summary$family_id[1],
  adduct_fam = adduct_fam
)
} # }
```
