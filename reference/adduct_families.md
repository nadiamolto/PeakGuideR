# Build neutral-mass-based adduct families from pairwise adduct candidate edges

Groups adduct candidate edges into putative adduct families using the
pairwise output of
[`adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_candidates.md).
In contrast to simple graph-based grouping, this function groups edges
by their inferred neutral mass, so that alternative adduct hypotheses
involving the same feature are kept as separate neutral-mass families
rather than being forced into the same connected component.

The function returns:

- a family-level summary table (`family_summary`)

- a feature-level membership table (`family_members`)

- an edge-level table with family assignment (`family_edges`)

This is intended to provide chemically more interpretable adduct
families, where each family represents a candidate neutral mass
supported by one or more mass-compatible and spatially correlated adduct
relationships.

## Usage

``` r
adduct_families(
  adduct_edges,
  min_score_adduct = 0.3,
  neutral_cluster_ppm = 5,
  min_family_size = 2L,
  use_only_valid = TRUE
)
```

## Arguments

- adduct_edges:

  A `data.frame` produced by
  [`adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_candidates.md).
  Must contain at least: `idx_i`, `idx_j`, `mz_i`, `mz_j`,
  `score_adduct`, `neutral_mass_i`, `neutral_mass_j`,
  `neutral_mass_mean`, `adduct_i`, `adduct_j`, `is_valid_adduct`.

- min_score_adduct:

  Numeric in the range 0 to 1. Minimum adduct score required to include
  an edge. Default `0.2`.

- neutral_cluster_ppm:

  Numeric. PPM tolerance used to cluster edges into neutral-mass
  families. Default `5`.

- min_family_size:

  Integer. Minimum number of unique features required to keep a family.
  Default `2L`.

- use_only_valid:

  Logical. If `TRUE`, only edges with `is_valid_adduct == TRUE` are
  used. Default `TRUE`.

## Value

A list with:

- `family_summary`: one row per neutral-mass adduct family

- `family_members`: one row per feature/adduct role within a family

- `family_edges`: one row per adduct edge assigned to a family

`family_summary` contains:

- `family_id`

- `neutral_mass_consensus`

- `neutral_mass_range_ppm`

- `family_size`

- `n_edges`

- `mean_score_adduct`

- `median_score_adduct`

- `adducts`

- `feature_idx`

- `has_role_conflict`

`family_members` contains:

- `family_id`

- `idx`

- `mz`

- `adduct`

- `neutral_mass_from_feature`

- `n_edges`

- `mean_edge_score`

`family_edges` contains the filtered adduct candidate edges plus:

- `family_id`

## Examples

``` r
if (FALSE) { # \dontrun{
adduct_res <- adduct_candidates(pkm, ion_mode = "pos")

adduct_fam <- adduct_families(
  adduct_res,
  min_score_adduct = 0.3,
  neutral_cluster_ppm = 5
)

head(adduct_fam$family_summary)
head(adduct_fam$family_members)
head(adduct_fam$family_edges)
} # }
```
