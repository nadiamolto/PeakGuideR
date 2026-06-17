#' Build neutral-mass-based adduct families from pairwise adduct candidate edges
#'
#' @description
#' Groups adduct candidate edges into putative adduct families using the pairwise
#' output of `adduct_candidates()`. In contrast to simple graph-based grouping,
#' this function groups edges by their inferred neutral mass, so that alternative
#' adduct hypotheses involving the same feature are kept as separate neutral-mass
#' families rather than being forced into the same connected component.
#'
#' The function returns:
#' \itemize{
#'   \item a family-level summary table (`family_summary`)
#'   \item a feature-level membership table (`family_members`)
#'   \item an edge-level table with family assignment (`family_edges`)
#' }
#'
#' This is intended to provide chemically more interpretable adduct families,
#' where each family represents a candidate neutral mass supported by one or more
#' mass-compatible and spatially correlated adduct relationships.
#'
#' @param adduct_edges A `data.frame` produced by `adduct_candidates()`.
#'   Must contain at least:
#'   `idx_i`, `idx_j`, `mz_i`, `mz_j`, `score_adduct`,
#'   `neutral_mass_i`, `neutral_mass_j`, `neutral_mass_mean`,
#'   `adduct_i`, `adduct_j`, `is_valid_adduct`.
#' @param min_score_adduct Numeric in the range 0 to 1. Minimum adduct score
#'   required to include an edge. Default `0.2`.
#' @param neutral_cluster_ppm Numeric. PPM tolerance used to cluster edges into
#'   neutral-mass families. Default `5`.
#' @param min_family_size Integer. Minimum number of unique features required to
#'   keep a family. Default `2L`.
#' @param use_only_valid Logical. If `TRUE`, only edges with
#'   `is_valid_adduct == TRUE` are used. Default `TRUE`.
#'
#' @return A list with:
#' \itemize{
#'   \item `family_summary`: one row per neutral-mass adduct family
#'   \item `family_members`: one row per feature/adduct role within a family
#'   \item `family_edges`: one row per adduct edge assigned to a family
#' }
#'
#' `family_summary` contains:
#' \itemize{
#'   \item `family_id`
#'   \item `neutral_mass_consensus`
#'   \item `neutral_mass_range_ppm`
#'   \item `family_size`
#'   \item `n_edges`
#'   \item `mean_score_adduct`
#'   \item `median_score_adduct`
#'   \item `adducts`
#'   \item `feature_idx`
#'   \item `has_role_conflict`
#' }
#'
#' `family_members` contains:
#' \itemize{
#'   \item `family_id`
#'   \item `idx`
#'   \item `mz`
#'   \item `adduct`
#'   \item `neutral_mass_from_feature`
#'   \item `n_edges`
#'   \item `mean_edge_score`
#' }
#'
#' `family_edges` contains the filtered adduct candidate edges plus:
#' \itemize{
#'   \item `family_id`
#' }
#'
#' @examples
#' \dontrun{
#' adduct_res <- adduct_candidates(pkm, ion_mode = "pos")
#'
#' adduct_fam <- adduct_families(
#'   adduct_res,
#'   min_score_adduct = 0.3,
#'   neutral_cluster_ppm = 5
#' )
#'
#' head(adduct_fam$family_summary)
#' head(adduct_fam$family_members)
#' head(adduct_fam$family_edges)
#' }
#' @export
adduct_families <- function(
    adduct_edges,
    min_score_adduct = 0.3,
    neutral_cluster_ppm = 5,
    min_family_size = 2L,
    use_only_valid = TRUE
) {
  stopifnot(is.data.frame(adduct_edges))

  required_cols <- c(
    "idx_i", "idx_j",
    "mz_i", "mz_j",
    "score_adduct",
    "neutral_mass_i", "neutral_mass_j", "neutral_mass_mean",
    "adduct_i", "adduct_j",
    "is_valid_adduct"
  )
  stopifnot(all(required_cols %in% names(adduct_edges)))

  min_family_size <- as.integer(min_family_size)
  if (!is.finite(min_family_size) || min_family_size < 1L) {
    min_family_size <- 2L
  }

  neutral_cluster_ppm <- as.numeric(neutral_cluster_ppm)
  if (!is.finite(neutral_cluster_ppm) || neutral_cluster_ppm <= 0) {
    neutral_cluster_ppm <- 5
  }

  # 1) Filter edges
  edges <- adduct_edges[
    is.finite(adduct_edges$score_adduct) &
      adduct_edges$score_adduct >= min_score_adduct &
      is.finite(adduct_edges$neutral_mass_mean),
    , drop = FALSE
  ]

  if (isTRUE(use_only_valid)) {
    edges <- edges[
      edges$is_valid_adduct %in% TRUE,
      , drop = FALSE
    ]
  }

  if (!nrow(edges)) {
    empty_summary <- data.frame(
      family_id = integer(),
      neutral_mass_consensus = numeric(),
      neutral_mass_range_ppm = numeric(),
      family_size = integer(),
      n_edges = integer(),
      mean_score_adduct = numeric(),
      median_score_adduct = numeric(),
      adducts = character(),
      feature_idx = character(),
      has_role_conflict = logical(),
      stringsAsFactors = FALSE
    )

    empty_members <- data.frame(
      family_id = integer(),
      idx = integer(),
      mz = numeric(),
      adduct = character(),
      neutral_mass_from_feature = numeric(),
      n_edges = integer(),
      mean_edge_score = numeric(),
      stringsAsFactors = FALSE
    )

    empty_edges <- edges
    empty_edges$family_id <- integer()

    return(list(
      family_summary = empty_summary,
      family_members = empty_members,
      family_edges = empty_edges
    ))
  }

  # 2) Cluster edges by inferred neutral mass
  edges <- edges[order(edges$neutral_mass_mean), , drop = FALSE]
  rownames(edges) <- NULL

  family_id <- integer(nrow(edges))
  current_family <- 1L
  family_id[1L] <- current_family
  current_ref <- edges$neutral_mass_mean[1L]

  if (nrow(edges) >= 2L) {
    for (i in 2L:nrow(edges)) {
      ppm_diff <- 1e6 * abs(edges$neutral_mass_mean[i] - current_ref) /
        pmax(abs(current_ref), 1e-12)

      if (is.finite(ppm_diff) && ppm_diff <= neutral_cluster_ppm) {
        family_id[i] <- current_family
        current_ref <- stats::median(
          edges$neutral_mass_mean[family_id == current_family],
          na.rm = TRUE
        )
      } else {
        current_family <- current_family + 1L
        family_id[i] <- current_family
        current_ref <- edges$neutral_mass_mean[i]
      }
    }
  }

  edges$family_id <- family_id

  # 3) Remove families smaller than min_family_size
  family_sizes <- vapply(
    split(edges, edges$family_id),
    function(x) length(unique(c(x$idx_i, x$idx_j))),
    integer(1)
  )

  keep_families <- as.integer(names(family_sizes)[family_sizes >= min_family_size])

  edges <- edges[
    edges$family_id %in% keep_families,
    , drop = FALSE
  ]

  if (!nrow(edges)) {
    empty_summary <- data.frame(
      family_id = integer(),
      neutral_mass_consensus = numeric(),
      neutral_mass_range_ppm = numeric(),
      family_size = integer(),
      n_edges = integer(),
      mean_score_adduct = numeric(),
      median_score_adduct = numeric(),
      adducts = character(),
      feature_idx = character(),
      has_role_conflict = logical(),
      stringsAsFactors = FALSE
    )

    empty_members <- data.frame(
      family_id = integer(),
      idx = integer(),
      mz = numeric(),
      adduct = character(),
      neutral_mass_from_feature = numeric(),
      n_edges = integer(),
      mean_edge_score = numeric(),
      stringsAsFactors = FALSE
    )

    return(list(
      family_summary = empty_summary,
      family_members = empty_members,
      family_edges = edges
    ))
  }

  # Re-number families consecutively after filtering
  old_ids <- sort(unique(edges$family_id))
  new_ids <- stats::setNames(seq_along(old_ids), old_ids)
  edges$family_id <- unname(new_ids[as.character(edges$family_id)])

  # 4) Build family-level and member-level outputs
  member_rows <- list()
  summary_rows <- list()

  for (g in sort(unique(edges$family_id))) {
    edge_sub <- edges[
      edges$family_id == g,
      , drop = FALSE
    ]

    neutral_pool <- edge_sub$neutral_mass_mean
    neutral_pool <- neutral_pool[is.finite(neutral_pool)]

    neutral_consensus <- if (length(neutral_pool)) {
      stats::median(neutral_pool, na.rm = TRUE)
    } else {
      NA_real_
    }

    neutral_range <- if (length(neutral_pool) >= 2L) {
      diff(range(neutral_pool, na.rm = TRUE))
    } else {
      0
    }

    neutral_range_ppm <- if (is.finite(neutral_consensus) && neutral_consensus != 0) {
      1e6 * neutral_range / abs(neutral_consensus)
    } else {
      NA_real_
    }

    # Build feature/adduct-role member table
    members_i <- data.frame(
      family_id = g,
      idx = edge_sub$idx_i,
      mz = edge_sub$mz_i,
      adduct = edge_sub$adduct_i,
      neutral_mass_from_feature = edge_sub$neutral_mass_i,
      score_adduct = edge_sub$score_adduct,
      stringsAsFactors = FALSE
    )

    members_j <- data.frame(
      family_id = g,
      idx = edge_sub$idx_j,
      mz = edge_sub$mz_j,
      adduct = edge_sub$adduct_j,
      neutral_mass_from_feature = edge_sub$neutral_mass_j,
      score_adduct = edge_sub$score_adduct,
      stringsAsFactors = FALSE
    )

    members_raw <- rbind(members_i, members_j)

    member_key <- paste(
      members_raw$family_id,
      members_raw$idx,
      members_raw$adduct,
      sep = "__"
    )

    member_split <- split(members_raw, member_key)

    members <- do.call(
      rbind,
      lapply(member_split, function(z) {
        data.frame(
          family_id = z$family_id[1L],
          idx = z$idx[1L],
          mz = stats::median(z$mz, na.rm = TRUE),
          adduct = z$adduct[1L],
          neutral_mass_from_feature = stats::median(
            z$neutral_mass_from_feature,
            na.rm = TRUE
          ),
          n_edges = nrow(z),
          mean_edge_score = mean(z$score_adduct, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )

    rownames(members) <- NULL

    # A role conflict occurs when the same feature appears with more than one
    # adduct role within the same neutral-mass family.
    role_count <- stats::aggregate(
      adduct ~ idx,
      data = members,
      FUN = function(x) length(unique(x))
    )

    has_role_conflict <- any(role_count$adduct > 1L)

    feature_idx <- paste(sort(unique(members$idx)), collapse = ",")
    adducts <- paste(sort(unique(members$adduct)), collapse = ",")

    member_rows[[as.character(g)]] <- members

    summary_rows[[as.character(g)]] <- data.frame(
      family_id = g,
      neutral_mass_consensus = neutral_consensus,
      neutral_mass_range_ppm = neutral_range_ppm,
      family_size = length(unique(members$idx)),
      n_edges = nrow(edge_sub),
      mean_score_adduct = mean(edge_sub$score_adduct, na.rm = TRUE),
      median_score_adduct = stats::median(edge_sub$score_adduct, na.rm = TRUE),
      adducts = adducts,
      feature_idx = feature_idx,
      has_role_conflict = has_role_conflict,
      stringsAsFactors = FALSE
    )
  }

  family_members <- do.call(rbind, member_rows)
  rownames(family_members) <- NULL

  family_summary <- do.call(rbind, summary_rows)
  rownames(family_summary) <- NULL

  family_members <- family_members[
    order(family_members$family_id, family_members$mz, family_members$adduct),
    , drop = FALSE
  ]

  family_summary <- family_summary[
    order(
      -family_summary$family_size,
      -family_summary$mean_score_adduct,
      family_summary$neutral_mass_consensus
    ),
    , drop = FALSE
  ]

  rownames(family_summary) <- NULL

  list(
    family_summary = family_summary,
    family_members = family_members,
    family_edges = edges
  )
}

# adduct_fam <- adduct_families(
#  adduct_res_pos,
#  min_score_adduct = 0.3,
#  neutral_cluster_ppm = 5,
#  min_family_size = 2,
#  use_only_valid = TRUE
#)
