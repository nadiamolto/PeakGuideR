#' Build single-adduct candidate hypotheses for unassigned features
#'
#' @description
#' Builds candidate neutral-mass hypotheses for m/z features that were not
#' assigned to any adduct family by [adduct_families()], excluding features
#' already explained as isotope satellite peaks (`is_c13_m1`, `is_c13_m2`,
#' `is_eips_isotope` in `feature_summary`).
#'
#' For each unassigned feature and each adduct compatible with `ion_mode`, a
#' neutral-mass hypothesis is inferred (`neutral_mass = mz - adduct_mass`,
#' the same convention as elsewhere in PeakGuideR). These hypotheses are then
#' clustered by ppm tolerance using the same mechanism as
#' [adduct_families()] / [build_neutral_mass_candidates()]
#' (see [assign_neutral_mass_clusters()]), and share the **same**
#' `neutral_mass_id` space as adduct families: a hypothesis that falls within
#' `neutral_cluster_ppm` of an already-detected family's neutral mass is
#' folded into that family's `neutral_mass_id` (`hypothesis_origin =
#' "single_adduct_clustered"`) instead of becoming a separate entry. Isolated
#' hypotheses that do not fall into any existing family are clustered only
#' among themselves and receive new `neutral_mass_id` values
#' (`hypothesis_origin = "single_adduct_isolated"`).
#'
#' An `adduct_spatial_score` is computed for every hypothesis by correlating
#' the feature's intensity image against every other feature already linked
#' to the same `neutral_mass_id` (family members and/or other single-adduct
#' hypotheses), reusing the spatial-similarity utilities in
#' `R/spatial_similarity.R`. A hypothesis with no such partner (a truly
#' isolated cluster of size one) gets `adduct_spatial_score = NA`, not `0`.
#'
#' Candidates are then retrieved by mass from `compound_db`
#' ([match_compound_db_by_mass()]) and, when `matrix = "HCCA"`, from
#' `standards_db` ([match_standards_by_mass()]), and identity is resolved
#' with [resolve_cross_source_identity()] exactly as for family-derived
#' candidates. `standard_adduct_recovery_score` is computed against the
#' detected-adducts set enriched by the neutral-mass clustering above (family
#' members' adducts and every single-adduct hypothesis sharing the same
#' `neutral_mass_id`).
#'
#' @param pkm A peak matrix list with `mass` and `intensity`.
#' @param feature_summary Output from [build_feature_summary()].
#' @param adduct_fam Output from [adduct_families()].
#' @param ion_mode `"pos"` or `"neg"`.
#' @param matrix Matrix name. Standard-adduct support is applied only when
#'   `matrix = "HCCA"`. Use `NULL` to skip it.
#' @param adducts Optional adduct definition data.frame. If `NULL`,
#'   `default_adducts(ion_mode)` is used.
#' @param compound_db Compound mass database. If `NULL`, the included example
#'   database is loaded.
#' @param standards_db Standard adduct library. If `NULL` and `matrix =
#'   "HCCA"`, the included example library is loaded.
#' @param ppm_tol PPM tolerance for compound/standard mass matching.
#' @param neutral_cluster_ppm PPM tolerance used to cluster inferred neutral
#'   masses, including folding hypotheses into existing adduct families.
#' @param top_n Maximum number of compound candidates kept per neutral mass.
#'   Use `NULL` to retain all candidates within `ppm_tol`.
#' @param method Spatial similarity metric: `"pearson"`, `"cosine"` or
#'   `"spearman"`.
#' @param transform Intensity transformation before similarity calculation:
#'   `"none"`, `"log1p"` or `"zscore"`.
#' @param min_quantile Numeric in the range 0 to 1. Feature-wise low-intensity
#'   quantile filter used by the spatial similarity score.
#' @param clip_negatives Logical. If `TRUE`, negative intensities are
#'   truncated to zero before transformation.
#' @param quiet Logical. If `FALSE`, database loading functions may print
#'   notices.
#'
#' @return A data.frame with one row per feature, inferred adduct and
#'   candidate identity: `neutral_mass_id`, `neutral_mass_consensus`,
#'   `feature_idx`, `inferred_adduct`, `hypothesis_origin`,
#'   `adduct_spatial_score`, `source`, `identity_match_type`,
#'   `possible_link_formula_only`, `candidate_ppm_error`,
#'   `candidate_neutral_mass`, `broad_db_*`, `standard_db_*` and
#'   `standard_adduct_recovery_score`.
#'
#' @examples
#' \dontrun{
#' data(example_pkm, package = "PeakGuideR")
#'
#' morph_results <- iso_morphology_candidates(example_pkm, prefer_mode = "ppm")
#' cir_results <- cir_score(morph_results, example_pkm)
#' eips_results <- eips_score(
#'   morph_results, example_pkm,
#'   ion_mode = "pos", cir_df = cir_results, morph_df = morph_results
#' )
#'
#' adduct_edges <- adduct_candidates(example_pkm, ion_mode = "pos")
#' adduct_fam <- adduct_families(adduct_edges)
#'
#' relation_table <- build_relation_table(
#'   cir_results = cir_results,
#'   eips_results = eips_results,
#'   adduct_fam = adduct_fam
#' )
#' feature_summary <- build_feature_summary(
#'   relation_table = relation_table,
#'   adduct_fam = adduct_fam,
#'   pkm = example_pkm
#' )
#'
#' single_adduct_hyp <- build_single_adduct_candidates(
#'   pkm = example_pkm,
#'   feature_summary = feature_summary,
#'   adduct_fam = adduct_fam,
#'   ion_mode = "pos",
#'   matrix = "HCCA"
#' )
#' head(single_adduct_hyp)
#' }
#' @export
build_single_adduct_candidates <- function(
    pkm,
    feature_summary,
    adduct_fam,
    ion_mode = c("pos", "neg"),
    matrix = NULL,
    adducts = NULL,
    compound_db = NULL,
    standards_db = NULL,
    ppm_tol = 5,
    neutral_cluster_ppm = 5,
    top_n = 10L,
    method = c("pearson", "cosine", "spearman"),
    transform = c("none", "log1p", "zscore"),
    min_quantile = 0.01,
    clip_negatives = TRUE,
    quiet = FALSE
) {
  ion_mode <- match.arg(ion_mode)
  method <- match.arg(method)
  transform <- match.arg(transform)

  if (!is.list(pkm) || is.null(pkm$mass) || is.null(pkm$intensity)) {
    stop("`pkm` must be a list containing `mass` and `intensity`.", call. = FALSE)
  }
  if (!is.matrix(pkm$intensity) || !is.numeric(pkm$intensity)) {
    stop("`pkm$intensity` must be a numeric matrix.", call. = FALSE)
  }
  if (length(pkm$mass) != ncol(pkm$intensity)) {
    stop(
      "Length of `pkm$mass` must match the number of columns in `pkm$intensity`.",
      call. = FALSE
    )
  }

  stopifnot(is.data.frame(feature_summary))
  if (!all(c("idx", "mz") %in% names(feature_summary))) {
    stop("`feature_summary` must contain `idx` and `mz`.", call. = FALSE)
  }

  stopifnot(is.list(adduct_fam))
  if (!all(c("family_summary", "family_members") %in% names(adduct_fam))) {
    stop(
      "`adduct_fam` must contain `family_summary` and `family_members`.",
      call. = FALSE
    )
  }

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1L) {
      stop("`top_n` must be NULL or a positive integer.", call. = FALSE)
    }
  }

  ppm_tol <- as.numeric(ppm_tol)
  if (!is.finite(ppm_tol) || ppm_tol <= 0) ppm_tol <- 5

  neutral_cluster_ppm <- as.numeric(neutral_cluster_ppm)
  if (!is.finite(neutral_cluster_ppm) || neutral_cluster_ppm <= 0) {
    neutral_cluster_ppm <- ppm_tol
  }

  empty_result <- function() {
    data.frame(
      neutral_mass_id = integer(),
      neutral_mass_consensus = numeric(),
      feature_idx = integer(),
      inferred_adduct = character(),
      hypothesis_origin = character(),
      adduct_spatial_score = numeric(),
      source = character(),
      identity_match_type = character(),
      possible_link_formula_only = character(),
      candidate_ppm_error = numeric(),
      candidate_neutral_mass = numeric(),
      broad_db_name = character(),
      broad_db_formula = character(),
      broad_db_inchikey = character(),
      broad_db_inchi = character(),
      broad_db_smiles = character(),
      broad_db_id = character(),
      broad_db_source = character(),
      standard_db_name = character(),
      standard_db_formula = character(),
      standard_db_inchikey = character(),
      standard_db_smiles = character(),
      standard_db_compound_id = character(),
      standard_db_adducts = character(),
      standard_adduct_recovery_score = numeric(),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(adducts)) adducts <- default_adducts(ion_mode)
  adf <- adducts[adducts$mode == ion_mode, , drop = FALSE]
  if (!nrow(adf)) {
    stop("No adducts available for ion_mode = '", ion_mode, "'.", call. = FALSE)
  }

  fs <- feature_summary
  for (col in c("has_adduct_family", "is_c13_m1", "is_c13_m2", "is_eips_isotope")) {
    if (!col %in% names(fs)) fs[[col]] <- FALSE
    fs[[col]] <- dplyr::coalesce(fs[[col]], FALSE)
  }

  unassigned <- fs |>
    dplyr::filter(!has_adduct_family, !is_c13_m1, !is_c13_m2, !is_eips_isotope) |>
    dplyr::select(idx, mz) |>
    dplyr::distinct()

  if (!nrow(unassigned)) {
    return(empty_result())
  }

  hyp <- tidyr::crossing(
    unassigned,
    adf |> dplyr::select(adduct = name, adduct_mass = mass)
  ) |>
    dplyr::mutate(neutral_mass_hypothesis = mz - adduct_mass) |>
    dplyr::filter(is.finite(neutral_mass_hypothesis))

  if (!nrow(hyp)) {
    return(empty_result())
  }

  # --- Shared neutral_mass_id space with adduct families ---------------
  fam_summary <- adduct_fam$family_summary
  fam_id_map <- data.frame(family_id = integer(), neutral_mass_id = integer())
  family_clusters <- data.frame(neutral_mass_id = integer(), neutral_mass_consensus = numeric())

  if (is.data.frame(fam_summary) && nrow(fam_summary) &&
      all(c("family_id", "neutral_mass_consensus") %in% names(fam_summary))) {
    fam_mass_all <- as.numeric(fam_summary$neutral_mass_consensus)
    keep_fam <- is.finite(fam_mass_all)

    if (any(keep_fam)) {
      fam_cluster_id <- assign_neutral_mass_clusters(fam_mass_all[keep_fam], neutral_cluster_ppm)

      fam_id_map <- data.frame(
        family_id = fam_summary$family_id[keep_fam],
        neutral_mass_id = fam_cluster_id
      )

      family_clusters <- data.frame(
        neutral_mass_id = fam_cluster_id,
        neutral_mass_consensus = fam_mass_all[keep_fam]
      ) |>
        dplyr::group_by(neutral_mass_id) |>
        dplyr::summarise(neutral_mass_consensus = mean(neutral_mass_consensus), .groups = "drop")
    }
  }

  nearest_family_id <- function(mass_vec) {
    if (!nrow(family_clusters)) {
      return(rep(NA_integer_, length(mass_vec)))
    }
    vapply(mass_vec, function(m) {
      ppm_diff <- 1e6 * abs(family_clusters$neutral_mass_consensus - m) /
        pmax(abs(m), 1e-12)
      best <- which.min(ppm_diff)
      if (length(best) == 1L && is.finite(ppm_diff[best]) &&
          ppm_diff[best] <= neutral_cluster_ppm) {
        family_clusters$neutral_mass_id[best]
      } else {
        NA_integer_
      }
    }, integer(1))
  }

  hyp$.family_cluster <- nearest_family_id(hyp$neutral_mass_hypothesis)

  isolated_mask <- is.na(hyp$.family_cluster)
  next_id <- if (nrow(family_clusters)) max(family_clusters$neutral_mass_id) else 0L

  hyp$.isolated_cluster <- NA_integer_
  if (any(isolated_mask)) {
    iso_cluster <- assign_neutral_mass_clusters(
      hyp$neutral_mass_hypothesis[isolated_mask],
      neutral_cluster_ppm
    )
    hyp$.isolated_cluster[isolated_mask] <- iso_cluster + next_id
  }

  hyp$neutral_mass_id <- dplyr::coalesce(hyp$.family_cluster, hyp$.isolated_cluster)
  hyp$hypothesis_origin <- dplyr::if_else(
    !is.na(hyp$.family_cluster),
    "single_adduct_clustered",
    "single_adduct_isolated"
  )

  iso_consensus <- hyp |>
    dplyr::filter(hypothesis_origin == "single_adduct_isolated") |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(.iso_consensus = mean(neutral_mass_hypothesis), .groups = "drop")

  hyp <- hyp |>
    dplyr::left_join(
      family_clusters |> dplyr::rename(.fam_consensus = neutral_mass_consensus),
      by = "neutral_mass_id"
    ) |>
    dplyr::left_join(iso_consensus, by = "neutral_mass_id") |>
    dplyr::mutate(
      neutral_mass_consensus = dplyr::coalesce(.fam_consensus, .iso_consensus)
    )

  # --- Adduct spatial score ---------------------------------------------
  partners_by_nmid <- list()

  if (nrow(fam_id_map) && is.data.frame(adduct_fam$family_members) &&
      nrow(adduct_fam$family_members)) {
    fam_members_nm <- adduct_fam$family_members |>
      dplyr::inner_join(fam_id_map, by = "family_id") |>
      dplyr::select(neutral_mass_id, idx) |>
      dplyr::distinct()

    for (nmid in unique(fam_members_nm$neutral_mass_id)) {
      partners_by_nmid[[as.character(nmid)]] <- unique(
        fam_members_nm$idx[fam_members_nm$neutral_mass_id == nmid]
      )
    }
  }

  hyp_idx_by_nmid <- split(hyp$idx, hyp$neutral_mass_id)
  for (nm in names(hyp_idx_by_nmid)) {
    partners_by_nmid[[nm]] <- unique(c(partners_by_nmid[[nm]], hyp_idx_by_nmid[[nm]]))
  }

  score_cache <- new.env(parent = emptyenv())

  pairwise_score <- function(i, j) {
    key <- paste(min(i, j), max(i, j), sep = "_")
    if (exists(key, envir = score_cache, inherits = FALSE)) {
      return(get(key, envir = score_cache, inherits = FALSE))
    }
    pp <- preprocess_xy(
      pkm$intensity[, i], pkm$intensity[, j],
      min_quantile = min_quantile,
      clip_negatives = clip_negatives,
      transform = transform
    )
    s <- if (is.null(pp)) NA_real_ else score_core(pp$x, pp$y, method = method)
    assign(key, s, envir = score_cache)
    s
  }

  hyp$adduct_spatial_score <- vapply(seq_len(nrow(hyp)), function(r) {
    nm <- as.character(hyp$neutral_mass_id[r])
    this_idx <- hyp$idx[r]
    partners <- setdiff(partners_by_nmid[[nm]], this_idx)
    if (!length(partners)) return(NA_real_)
    scores <- vapply(partners, function(p) pairwise_score(this_idx, p), numeric(1))
    if (all(is.na(scores))) return(NA_real_)
    max(scores, na.rm = TRUE)
  }, numeric(1))

  # --- Candidate retrieval + identity resolution -------------------------
  neutral_for_matching <- hyp |>
    dplyr::select(neutral_mass_id, neutral_mass_consensus) |>
    dplyr::distinct()

  if (is.null(compound_db)) {
    compound_db <- load_compound_mass_database(quiet = quiet)
  }
  stopifnot(is.data.frame(compound_db))

  compound_candidates <- match_compound_db_by_mass(
    neutral_for_matching, compound_db,
    ppm_tol = ppm_tol, top_n = top_n
  )

  standard_candidates <- NULL
  if (!is.null(matrix) && identical(toupper(matrix), "HCCA")) {
    if (is.null(standards_db)) {
      standards_db <- load_standards_adduct_library(quiet = quiet)
    }
    stopifnot(is.data.frame(standards_db))

    standard_candidates <- match_standards_by_mass(
      neutral_for_matching, standards_db,
      ion_mode = ion_mode, matrix = matrix, ppm_tol = ppm_tol
    )
  }

  if (!nrow(compound_candidates) &&
      (is.null(standard_candidates) || !nrow(standard_candidates))) {
    resolved_identity <- data.frame(
      neutral_mass_id = integer(),
      source = character(),
      identity_match_type = character(),
      possible_link_formula_only = character(),
      candidate_ppm_error = numeric(),
      candidate_neutral_mass = numeric(),
      broad_db_name = character(),
      broad_db_formula = character(),
      broad_db_inchikey = character(),
      broad_db_inchi = character(),
      broad_db_smiles = character(),
      broad_db_id = character(),
      broad_db_source = character(),
      standard_db_name = character(),
      standard_db_formula = character(),
      standard_db_inchikey = character(),
      standard_db_smiles = character(),
      standard_db_compound_id = character(),
      standard_db_adducts = character(),
      stringsAsFactors = FALSE
    )
  } else {
    resolved_identity <- resolve_cross_source_identity(compound_candidates, standard_candidates)
  }

  # --- Detected-adducts set enriched by the shared clustering -------------
  family_adducts_nm <- NULL
  if (nrow(fam_id_map) && is.data.frame(adduct_fam$family_members) &&
      nrow(adduct_fam$family_members)) {
    family_adducts_nm <- adduct_fam$family_members |>
      dplyr::inner_join(fam_id_map, by = "family_id") |>
      dplyr::select(neutral_mass_id, adduct)
  }

  detected_map <- compute_detected_adducts_map(
    family_adducts_nm,
    hyp |> dplyr::select(neutral_mass_id, adduct)
  )

  out <- hyp |>
    dplyr::transmute(
      neutral_mass_id,
      neutral_mass_consensus,
      feature_idx = idx,
      inferred_adduct = adduct,
      hypothesis_origin,
      adduct_spatial_score
    ) |>
    dplyr::left_join(
      resolved_identity,
      by = "neutral_mass_id",
      relationship = "many-to-many"
    )

  out$standard_adduct_recovery_score <- mapply(
    function(nmid, std_adducts) {
      if (is.na(std_adducts)) return(NA_real_)
      det <- detected_map[[as.character(nmid)]]
      if (is.null(det)) det <- character(0)
      compute_standard_adduct_recovery(std_adducts, det)
    },
    out$neutral_mass_id, out$standard_db_adducts
  )

  out |>
    dplyr::arrange(neutral_mass_id, feature_idx, inferred_adduct, candidate_ppm_error)
}
