# Plot a peak evidence network

Plots a simple network of peak-to-peak relationships from a relation
table.

Nodes are MSI features and edges represent isotope, EIPS or adduct
evidence. The function does not perform annotation decisions; it only
visualizes the relations already present in `relation_table`.

## Usage

``` r
plot_peak_network(
  relation_table,
  min_score = 0.3,
  idx_focus = NULL,
  only_valid = TRUE,
  show_edge_labels = TRUE
)
```

## Arguments

- relation_table:

  Output from
  [`build_relation_table()`](https://nadiamolto.github.io/PeakGuideR/reference/build_relation_table.md).

- min_score:

  Minimum evidence score required to plot an edge.

- idx_focus:

  Optional vector of feature indices to focus on. If supplied, only
  relations involving these features and their direct neighbours are
  shown.

- only_valid:

  Logical. If `TRUE`, only valid relations are plotted.

- show_edge_labels:

  Logical. If `TRUE`, shows edge labels.

## Value

Invisibly returns the igraph object.
