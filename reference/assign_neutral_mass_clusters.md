# Cluster neutral masses by ppm tolerance

Shared internal sequential clustering used to assign `neutral_mass_id`
groups from a numeric vector of inferred neutral masses. Masses are
processed in ascending order; a new cluster starts whenever the ppm
distance to the current cluster's anchor mass (its first, lowest member)
exceeds `neutral_cluster_ppm`.

Used by both
[`build_neutral_mass_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_neutral_mass_candidates.md)
and
[`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)
so that, given the same adduct-family input and `neutral_cluster_ppm`,
they assign identical `neutral_mass_id` values to family-derived neutral
masses.

## Usage

``` r
assign_neutral_mass_clusters(masses, neutral_cluster_ppm)
```

## Arguments

- masses:

  Numeric vector of neutral masses.

- neutral_cluster_ppm:

  PPM tolerance.

## Value

Integer vector of cluster ids, same length and order as `masses`.
