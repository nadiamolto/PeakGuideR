# Plot a peak evidence network

Creates an interactive network of peak-to-peak relationships using
`visNetwork`.

...

Creates an interactive network of peak-to-peak relationships using
`visNetwork`.

Nodes are MSI features and edges represent isotope, EIPS or adduct
evidence. If a `peakguider_workflow` object is supplied, adduct-family
information is used when available to display adduct roles, family IDs
and consensus neutral masses.

## Usage

``` r
plot_peak_network(
  relation_table,
  min_score = 0.3,
  idx_focus = NULL,
  family_id = NULL,
  only_valid = TRUE,
  show_edge_labels = TRUE
)

plot_peak_network(
  relation_table,
  min_score = 0.3,
  idx_focus = NULL,
  family_id = NULL,
  only_valid = TRUE,
  show_edge_labels = TRUE
)
```

## Arguments

- relation_table:

  A `peakguider_workflow` object, a relation table from
  [`build_relation_table()`](https://nadiamolto.github.io/PeakGuideR/reference/build_relation_table.md),
  or an adduct edge table.

- min_score:

  Minimum evidence score required to plot an edge.

- idx_focus:

  Optional vector of feature indices to focus on. If supplied, only
  relations involving these features and their direct neighbours are
  shown.

- family_id:

  Optional adduct family ID to display. If supplied, the network is
  restricted to the selected adduct family when family information is
  available.

- only_valid:

  Logical. If `TRUE`, only valid relations are plotted when a validity
  column is available.

- show_edge_labels:

  Logical. If `TRUE`, shows edge labels.

## Value

A `visNetwork` htmlwidget.
