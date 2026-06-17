#' Plot a peak evidence network
#'
#' @description
#' Plots a simple network of peak-to-peak relationships from a relation table.
#'
#' Nodes are MSI features and edges represent isotope, EIPS or adduct evidence.
#' The function does not perform annotation decisions; it only visualizes the
#' relations already present in `relation_table`.
#'
#' @param relation_table Output from `build_relation_table()`.
#' @param min_score Minimum evidence score required to plot an edge.
#' @param idx_focus Optional vector of feature indices to focus on. If supplied,
#'   only relations involving these features and their direct neighbours are shown.
#' @param only_valid Logical. If `TRUE`, only valid relations are plotted.
#' @param show_edge_labels Logical. If `TRUE`, shows edge labels.
#'
#' @return Invisibly returns the igraph object.
#' @export
plot_peak_network <- function(
    relation_table,
    min_score = 0.3,
    idx_focus = NULL,
    only_valid = TRUE,
    show_edge_labels = TRUE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required to use `plot_peak_network()`. ",
      "Please install it with install.packages('igraph').",
      call. = FALSE
    )
  }

  stopifnot(is.data.frame(relation_table))

  if (!nrow(relation_table)) {
    message("Empty relation table.")
    return(invisible(NULL))
  }

  edges <- relation_table |>
    dplyr::filter(
      is.finite(evidence_score),
      evidence_score >= min_score
    )

  if (isTRUE(only_valid) && "is_valid" %in% names(edges)) {
    edges <- edges |>
      dplyr::filter(is_valid %in% TRUE)
  }

  if (!is.null(idx_focus)) {
    idx_focus <- as.integer(idx_focus)

    keep_idx <- unique(c(
      idx_focus,
      edges$from_idx[edges$to_idx %in% idx_focus],
      edges$to_idx[edges$from_idx %in% idx_focus]
    ))

    edges <- edges |>
      dplyr::filter(
        from_idx %in% keep_idx,
        to_idx %in% keep_idx
      )
  }

  if (!nrow(edges)) {
    message("No relations to plot after filtering.")
    return(invisible(NULL))
  }

  graph_edges <- edges |>
    dplyr::transmute(
      from = as.character(from_idx),
      to = as.character(to_idx),
      evidence_type = evidence_type,
      relation_type = relation_type,
      evidence_score = evidence_score,
      label = relation_type
    )

  g <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = FALSE
  )

  edge_labels <- if (isTRUE(show_edge_labels)) {
    igraph::E(g)$label
  } else {
    NA
  }

  graphics::plot(
    g,
    vertex.label = igraph::V(g)$name,
    vertex.size = 22,
    edge.width = 1 + 4 * igraph::E(g)$evidence_score,
    edge.label = edge_labels,
    main = "PeakGuideR evidence network"
  )

  invisible(g)
}
