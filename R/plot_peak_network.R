utils::globalVariables(c(
  ".data",
  ".score",
  "neutral_mass_mean",
  "neutral_mass_from_feature",
  "id",
  "label",
  "group",
  "from",
  "to",
  "title"
))

#' Plot a peak evidence network
#'
#' @description
#' Creates an interactive network of peak-to-peak relationships using
#' `visNetwork`.
#'
#' ...
plot_peak_network <- function(
    ...
) {
  ...
}
#' Plot a peak evidence network
#'
#' @description
#' Creates an interactive network of peak-to-peak relationships using
#' `visNetwork`.
#'
#' Nodes are MSI features and edges represent isotope, EIPS or adduct evidence.
#' If a `peakguider_workflow` object is supplied, adduct-family information is
#' used when available to display adduct roles, family IDs and consensus neutral
#' masses.
#'
#' @param relation_table A `peakguider_workflow` object, a relation table from
#'   `build_relation_table()`, or an adduct edge table.
#' @param min_score Minimum evidence score required to plot an edge.
#' @param idx_focus Optional vector of feature indices to focus on. If supplied,
#'   only relations involving these features and their direct neighbours are shown.
#' @param only_valid Logical. If `TRUE`, only valid relations are plotted when a
#'   validity column is available.
#' @param show_edge_labels Logical. If `TRUE`, shows edge labels.
#' @param family_id Optional adduct family ID to display. If supplied, the
#'   network is restricted to the selected adduct family when family information
#'   is available.
#'
#' @return A `visNetwork` htmlwidget.
#'
#' @export
plot_peak_network <- function(
    relation_table,
    min_score = 0.3,
    idx_focus = NULL,
    family_id = NULL,
    only_valid = TRUE,
    show_edge_labels = TRUE
) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop(
      "Package 'visNetwork' is required to use `plot_peak_network()`. ",
      "Please install it with install.packages('visNetwork').",
      call. = FALSE
    )
  }

  input <- relation_table

  family_members <- NULL
  family_summary <- NULL

  if (inherits(input, "peakguider_workflow")) {
    if (!is.null(input$adduct_families) &&
        is.list(input$adduct_families) &&
        "family_edges" %in% names(input$adduct_families) &&
        is.data.frame(input$adduct_families$family_edges) &&
        nrow(input$adduct_families$family_edges) > 0) {

      edges_raw <- input$adduct_families$family_edges
      edge_source <- "adduct_family_edges"

      family_members <- input$adduct_families$family_members
      family_summary <- input$adduct_families$family_summary

    } else if ("adduct_edges" %in% names(input) &&
               is.data.frame(input$adduct_edges) &&
               nrow(input$adduct_edges) > 0) {

      edges_raw <- input$adduct_edges
      edge_source <- "adduct_edges"

    } else if ("relation_table" %in% names(input) &&
               is.data.frame(input$relation_table)) {

      edges_raw <- input$relation_table
      edge_source <- "relation_table"

    } else {
      stop("No suitable edge or relation table found in `relation_table`.", call. = FALSE)
    }

  } else if (is.data.frame(input)) {
    edges_raw <- input

    if (all(c("idx_i", "idx_j") %in% names(edges_raw))) {
      edge_source <- "adduct_edges"
    } else {
      edge_source <- "relation_table"
    }

  } else {
    stop(
      "`relation_table` must be a peakguider_workflow object, a relation table, or an adduct edge table.",
      call. = FALSE
    )
  }

  if (!nrow(edges_raw)) {
    message("Empty relation table.")
    return(invisible(NULL))
  }

  # Edge table standardized

  if (all(c("idx_i", "idx_j") %in% names(edges_raw))) {
    # adduct_edges / family_edges style
    score_col <- dplyr::case_when(
      "score_adduct" %in% names(edges_raw) ~ "score_adduct",
      "score_spatial" %in% names(edges_raw) ~ "score_spatial",
      TRUE ~ NA_character_)

    if (is.na(score_col)) {
      edges <- edges_raw |>
        dplyr::mutate(.score = 1)
    } else {
      edges <- edges_raw |>
        dplyr::mutate(.score = .data[[score_col]])}

    edges <- edges |>
      dplyr::filter(is.finite(.score), .score >= min_score)

    if (isTRUE(only_valid) && "is_valid_adduct" %in% names(edges)) {
      edges <- edges |>
        dplyr::filter(is_valid_adduct %in% TRUE)
    }
    if (!is.null(family_id)) {
      family_id_filter <- as.integer(family_id)

      if ("family_id" %in% names(edges)) {
        edges <- edges |>
          dplyr::filter(.data$family_id %in% family_id_filter)

      } else if (!is.null(family_members) && is.data.frame(family_members)) {
        keep_family_idx <- family_members |>
          dplyr::filter(.data$family_id %in% family_id_filter) |>
          dplyr::pull(.data$idx) |>
          unique()

        edges <- edges |>
          dplyr::filter(
            .data$idx_i %in% keep_family_idx,
            .data$idx_j %in% keep_family_idx
          )
      }
    }
    if (!is.null(idx_focus)) {
      idx_focus <- as.integer(idx_focus)

      keep_idx <- unique(c(
        idx_focus,
        edges$idx_i[edges$idx_j %in% idx_focus],
        edges$idx_j[edges$idx_i %in% idx_focus]
      ))

      edges <- edges |>
        dplyr::filter(
          idx_i %in% keep_idx,
          idx_j %in% keep_idx)}

    if (!nrow(edges)) { message("No relations to plot after filtering.")
      return(invisible(NULL))}

    if (!"neutral_mass_mean" %in% names(edges)) {
      edges$neutral_mass_mean <- NA_real_
    }

    if (!"neutral_err_ppm" %in% names(edges)) {
      edges$neutral_err_ppm <- NA_real_
    }

    graph_edges <- edges |>
      dplyr::transmute(
        from = as.character(idx_i),
        to = as.character(idx_j),
        evidence_type = "adduct",
        relation_type = paste0(adduct_i, " - ", adduct_j),
        evidence_score = .score,
        neutral_mass_mean = as.numeric(neutral_mass_mean),
        label = relation_type,
        title = paste0(
          "<b>Adduct edge</b><br>",
          "Feature ", idx_i, " (", adduct_i, ", m/z ", signif(mz_i, 6), ")<br>",
          "Feature ", idx_j, " (", adduct_j, ", m/z ", signif(mz_j, 6), ")<br>",
          "Score: ", signif(.score, 3), "<br>",
          "Neutral mass: ", signif(neutral_mass_mean, 7), "<br>",
          "Neutral error ppm: ", signif(neutral_err_ppm, 3)
        )
      )

  } else {
    # relation_table style
    required_cols <- c("from_idx", "to_idx", "evidence_type", "relation_type", "evidence_score")
    missing_cols <- setdiff(required_cols, names(edges_raw))

    if (length(missing_cols) > 0) {
      stop(
        "Relation table is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    edges <- edges_raw |>
      dplyr::filter(
        is.finite(evidence_score),
        evidence_score >= min_score
      )

    if (isTRUE(only_valid) && "is_valid" %in% names(edges)) {
      edges <- edges |>
        dplyr::filter(is_valid %in% TRUE)
    }

    if (!is.null(family_id) && "group_id" %in% names(edges)) {
      family_id_filter <- as.integer(family_id)

      edges <- edges |>
        dplyr::filter(.data$group_id %in% family_id_filter)
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
      return(invisible(NULL))}

    graph_edges <- edges |>
      dplyr::transmute(
        from = as.character(from_idx),
        to = as.character(to_idx),
        evidence_type = evidence_type,
        relation_type = relation_type,
        evidence_score = evidence_score,
        neutral_mass_mean = NA_real_,
        label = relation_type,
        title = paste0(
          "<b>", evidence_type, "</b><br>",
          "Relation: ", relation_type, "<br>",
          "Feature ", from_idx,
          if ("from_mz" %in% names(edges)) paste0(" (m/z ", signif(from_mz, 6), ")") else "",
          "<br>",
          "Feature ", to_idx,
          if ("to_mz" %in% names(edges)) paste0(" (m/z ", signif(to_mz, 6), ")") else "",
          "<br>",
          "Score: ", signif(evidence_score, 3)
        )
      )
  }

  # Building noded

  node_ids <- sort(unique(c(graph_edges$from, graph_edges$to)))
  nodes <- data.frame(
    id = node_ids,
    label = node_ids,
    title = paste0("Feature ", node_ids),
    group = "feature",
    stringsAsFactors = FALSE
  )

  # Add family/adduct metadata (if it is available...)
  if (!is.null(family_members) && is.data.frame(family_members)) {
    fm <- family_members |>
      dplyr::filter(as.character(idx) %in% node_ids) |>
      dplyr::group_by(idx) |>
      dplyr::slice_max(order_by = mean_edge_score, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    if (nrow(fm) > 0) {
      nodes <- nodes |>
        dplyr::left_join(
          fm |>
            dplyr::transmute(
              id = as.character(idx),
              mz = mz,
              family_id = family_id,
              adduct = adduct,
              neutral_mass_from_feature = neutral_mass_from_feature,
              mean_edge_score = mean_edge_score
            ),
          by = "id"
        ) |>
        dplyr::mutate(
          label = dplyr::if_else(
            !is.na(adduct),
            paste0(id, "\n", adduct),
            label
          ),
          group = dplyr::if_else(
            !is.na(family_id),
            paste0("family_", family_id),
            group
          ),
          title = paste0(
            "<b>Feature ", id, "</b><br>",
            ifelse(!is.na(mz), paste0("m/z: ", signif(mz, 6), "<br>"), ""),
            ifelse(!is.na(family_id), paste0("Family: ", family_id, "<br>"), ""),
            ifelse(!is.na(adduct), paste0("Adduct: ", adduct, "<br>"), ""),
            ifelse(
              !is.na(neutral_mass_from_feature),
              paste0("Neutral mass from feature: ", signif(neutral_mass_from_feature, 7), "<br>"),
              ""
            ),
            ifelse(
              !is.na(mean_edge_score),
              paste0("Mean edge score: ", signif(mean_edge_score, 3)),
              ""
            )
          ))
    }}

  # Building visNetwork edges

  vis_edges <- graph_edges |>
    dplyr::transmute(
      from = from,
      to = to,
      label = if (isTRUE(show_edge_labels)) label else "",
      title = title,
      value = pmax(evidence_score, 0.05),
      width = 1 + 5 * pmax(evidence_score, 0),
      arrows = "")

  # ... and plotting them:

  visNetwork::visNetwork(nodes, vis_edges, height = "650px", width = "100%") |>
    visNetwork::visGroups(
      groupname = "feature",
      color = list(background = "#DDEBF7", border = "#3E6C99")
    ) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE,
      selectedBy = "group"
    ) |>
    visNetwork::visNodes(
      shape = "dot",
      size = 24,
      font = list(size = 18)
    ) |>
    visNetwork::visEdges(
      smooth = TRUE,
      font = list(size = 14, align = "middle")
    ) |>
    visNetwork::visLayout(
      randomSeed = 123
    ) |>
    visNetwork::visPhysics(
      solver = "forceAtlas2Based",
      stabilization = TRUE
    )
}
