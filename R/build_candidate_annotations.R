#' Default weights for the candidate priority score
#'
#' @description
#' Returns the default weights used by [compute_priority_score()] /
#' [build_candidate_annotations()] to combine per-candidate evidence scores
#' into a single `priority_score`.
#'
#' @return A named numeric vector of weights summing to 1:
#'   `mass_error_score` (0.20), `adduct_spatial_score` (0.20),
#'   `isotope_evidence_score` (0.20), `eips_evidence_score` (0.15),
#'   `standard_adduct_recovery_score` (0.15), `family_coherence_score`
#'   (0.10).
#'
#' @examples
#' default_priority_weights()
#'
#' @export
default_priority_weights <- function() {
  c(
    mass_error_score = 0.20,
    adduct_spatial_score = 0.20,
    isotope_evidence_score = 0.20,
    eips_evidence_score = 0.15,
    standard_adduct_recovery_score = 0.15,
    family_coherence_score = 0.10
  )
}


#' Compute a weighted candidate priority score
#'
#' @description
#' Combines per-candidate evidence scores into a single `priority_score` as a
#' weighted mean, excluding `NA` components from the calculation and
#' renormalizing the remaining weights to sum to 1. `NA` components are never
#' treated as `0`.
#'
#' @param scores A data.frame with one row per candidate and columns named
#'   after `weights` (missing columns are treated as all-`NA`), or a single
#'   named numeric vector for one candidate.
#' @param weights Named numeric vector of weights, see
#'   [default_priority_weights()].
#'
#' @return A numeric vector (or scalar) of priority scores in the range 0 to
#'   1, or `NA` where every component is `NA`.
#'
#' @examples
#' scores <- data.frame(
#'   mass_error_score = c(0.95, 0.80),
#'   adduct_spatial_score = c(0.70, NA),
#'   isotope_evidence_score = c(NA, NA),
#'   eips_evidence_score = c(NA, NA),
#'   standard_adduct_recovery_score = c(1, NA),
#'   family_coherence_score = c(0.65, 0.65)
#' )
#'
#' compute_priority_score(scores)
#' compute_priority_score(scores, weights = default_priority_weights())
#'
#' @export
compute_priority_score <- function(scores, weights = default_priority_weights()) {
  score_cols <- names(weights)

  if (is.data.frame(scores)) {
    mat_df <- as.data.frame(scores)
    for (col in score_cols) {
      if (!col %in% names(mat_df)) mat_df[[col]] <- NA_real_
    }
    mat <- as.matrix(mat_df[, score_cols, drop = FALSE])
  } else {
    scores <- as.list(scores)
    row <- vapply(score_cols, function(col) {
      if (is.null(scores[[col]])) NA_real_ else as.numeric(scores[[col]])
    }, numeric(1))
    mat <- matrix(row, nrow = 1, dimnames = list(NULL, score_cols))
  }

  storage.mode(mat) <- "double"
  w <- as.numeric(weights[score_cols])

  valid <- !is.na(mat)
  w_mat <- matrix(w, nrow = nrow(mat), ncol = length(w), byrow = TRUE)
  w_mat[!valid] <- 0

  mat_filled <- mat
  mat_filled[!valid] <- 0

  w_sum <- rowSums(w_mat)
  score_sum <- rowSums(mat_filled * w_mat)

  out <- ifelse(w_sum > 0, score_sum / w_sum, NA_real_)
  if (length(out) == 1L) unname(out) else unname(out)
}


#' Build the unified long-format candidate annotation table
#'
#' @description
#' Builds one long candidate-annotation table combining adduct-family
#' evidence and single-adduct hypotheses (see
#' [build_single_adduct_candidates()]), with database and standard-adduct
#' identity resolved via [resolve_cross_source_identity()] instead of the
#' identifier-hierarchy-free join used previously. This is the table returned
#' as `candidate_annotations` by [run_peakguider_workflow()].
#' [build_neutral_mass_candidates()] is a per-`neutral_mass_id` summary view
#' derived from this table.
#'
#' Each row represents one neutral mass, candidate identity (from
#' `compound_db`, `standards_db`, or both fused via a strong identifier
#' match) and evidence source. A match on molecular formula alone is never
#' fused into a single identity; see [resolve_cross_source_identity()].
#'
#' @param adduct_fam Output from [adduct_families()].
#' @param feature_summary Output from [build_feature_summary()].
#' @param pkm Optional peak matrix (`mass`, `intensity`). Required to include
#'   single-adduct hypotheses (see `include_single_adduct`); if `NULL`, the
#'   table contains family-derived candidates only.
#' @param ion_mode `"pos"` or `"neg"`.
#' @param matrix Matrix name. Standard-adduct support is applied only when
#'   `matrix = "HCCA"`. Use `NULL` to skip it.
#' @param adducts Optional adduct definition data.frame, used only for
#'   single-adduct hypotheses. If `NULL`, `default_adducts(ion_mode)` is
#'   used.
#' @param compound_db Compound mass database. If `NULL`, the included example
#'   database is loaded.
#' @param standards_db Standard adduct library. If `NULL` and `matrix =
#'   "HCCA"`, the included example library is loaded.
#' @param ppm_tol PPM tolerance for compound/standard mass matching.
#' @param neutral_cluster_ppm PPM tolerance used to cluster inferred neutral
#'   masses.
#' @param top_n Maximum number of compound candidates kept per neutral mass.
#'   Use `NULL` to retain all candidates within `ppm_tol`.
#' @param weights Named numeric vector of weights for `priority_score`; see
#'   [default_priority_weights()].
#' @param recovery_threshold Minimum `standard_adduct_recovery_score` (`NA`
#'   treated as `0`) required for `confidence_class =
#'   "cross_db_linked_adduct_recovered"` when `source == "both"`.
#' @param adduct_min_score Minimum `adduct_spatial_score` (`NA` treated as
#'   `0`) counted as multi-evidence support for `broad_db_only` candidates.
#'   Reuses the same threshold as the `adduct_min_score` workflow parameter.
#' @param ambiguity_gap Minimum `priority_score` gap between the top-1 and
#'   top-2 candidate of the same `neutral_mass_id` required to *not* flag
#'   `ambiguous_isomeric`.
#' @param include_single_adduct Logical. If `FALSE`, or if `pkm` is `NULL`,
#'   single-adduct hypotheses are not computed.
#' @param quiet Logical. If `FALSE`, database loading functions may print
#'   notices.
#'
#' @return A data.frame with one row per neutral mass, feature(s) and
#'   candidate identity, including `neutral_mass_id`,
#'   `neutral_mass_consensus`, `feature_idx`, `inferred_adduct`, `source`,
#'   `identity_match_type`, `possible_link_formula_only`,
#'   `hypothesis_origin`, `broad_db_*`, `standard_db_*`,
#'   `mass_error_score`, `adduct_spatial_score`, `isotope_evidence_score`,
#'   `eips_evidence_score`, `standard_adduct_recovery_score`,
#'   `family_coherence_score`, `priority_score`, `confidence_class` and
#'   `ambiguous_isomeric`.
#'
#' `confidence_class` summarizes the evidence combination behind each
#' candidate into one of six priority-ordered levels:
#' \itemize{
#'   \item `"cross_db_linked_adduct_recovered"`: `source == "both"` (the
#'     broad-database and standard-library identity were linked through a
#'     strong cross-database identifier) and `standard_adduct_recovery_score`
#'     (`NA` treated as `0`) is at least `recovery_threshold`. This reflects a
#'     cross-database identifier link plus recovered adduct evidence, not an
#'     analytical confirmation of the compound in the sample.
#'   \item `"cross_db_linked_partial_recovery"`: `source == "both"`, but the
#'     `recovery_threshold` above is not met.
#'   \item `"standards_only"`: `source == "standards_only"`.
#'   \item `"broad_db_only_multi_evidence"`: `source == "broad_db_only"` and
#'     at least two of `is_c13_m0`, `has_eips`,
#'     `adduct_spatial_score >= adduct_min_score` (`NA` treated as `0`) hold.
#'   \item `"broad_db_only_single_evidence"`: `source == "broad_db_only"` and
#'     exactly one of those three signals holds.
#'   \item `"broad_db_only_mass_only"`: `source == "broad_db_only"` and none
#'     of those three signals hold.
#' }
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
#' candidate_annotations <- build_candidate_annotations(
#'   adduct_fam = adduct_fam,
#'   feature_summary = feature_summary,
#'   pkm = example_pkm,
#'   ion_mode = "pos",
#'   matrix = "HCCA"
#' )
#' head(candidate_annotations)
#' }
#' @export
build_candidate_annotations <- function(
    adduct_fam,
    feature_summary,
    pkm = NULL,
    ion_mode = c("pos", "neg"),
    matrix = NULL,
    adducts = NULL,
    compound_db = NULL,
    standards_db = NULL,
    ppm_tol = 5,
    neutral_cluster_ppm = 5,
    top_n = 10L,
    weights = default_priority_weights(),
    recovery_threshold = 0.5,
    adduct_min_score = 0.5,
    ambiguity_gap = 0.05,
    include_single_adduct = TRUE,
    quiet = FALSE
) {
  ion_mode <- match.arg(ion_mode)

  stopifnot(is.list(adduct_fam))
  stopifnot(is.data.frame(feature_summary))

  if (!all(c("family_summary", "family_members") %in% names(adduct_fam))) {
    stop(
      "adduct_fam must contain `family_summary` and `family_members`.",
      call. = FALSE
    )
  }

  stopifnot(is.data.frame(adduct_fam$family_summary))
  stopifnot(is.data.frame(adduct_fam$family_members))

  required_family_cols <- c("family_id", "neutral_mass_consensus")
  missing_family_cols <- setdiff(required_family_cols, names(adduct_fam$family_summary))
  if (length(missing_family_cols) > 0) {
    stop(
      "adduct_fam$family_summary is missing required columns: ",
      paste(missing_family_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_member_cols <- c("family_id", "idx", "mz", "adduct")
  missing_member_cols <- setdiff(required_member_cols, names(adduct_fam$family_members))
  if (length(missing_member_cols) > 0) {
    stop(
      "adduct_fam$family_members is missing required columns: ",
      paste(missing_member_cols, collapse = ", "),
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

  if (is.null(compound_db)) {
    compound_db <- load_compound_mass_database(quiet = quiet)
  }
  stopifnot(is.data.frame(compound_db))
  if (!"MonoisotopicMass" %in% names(compound_db)) {
    stop("compound_db must contain a `MonoisotopicMass` column.", call. = FALSE)
  }

  use_standards <- !is.null(matrix) && identical(toupper(matrix), "HCCA")
  if (use_standards && is.null(standards_db)) {
    standards_db <- load_standards_adduct_library(quiet = quiet)
  }

  empty_annotations <- function() {
    data.frame(
      neutral_mass_id = integer(),
      neutral_mass_consensus = numeric(),
      feature_idx = character(),
      inferred_adduct = character(),
      source = character(),
      identity_match_type = character(),
      possible_link_formula_only = character(),
      hypothesis_origin = character(),
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
      mass_error_score = numeric(),
      adduct_spatial_score = numeric(),
      isotope_evidence_score = numeric(),
      eips_evidence_score = numeric(),
      standard_adduct_recovery_score = numeric(),
      family_coherence_score = numeric(),
      priority_score = numeric(),
      confidence_class = character(),
      ambiguous_isomeric = logical(),
      stringsAsFactors = FALSE
    )
  }

  # --- Family-level neutral_mass_id clustering (shared with
  # build_single_adduct_candidates()) --------------------------------------
  fam_summary <- adduct_fam$family_summary
  fam_members <- adduct_fam$family_members

  for (col in c("mean_score_adduct", "median_score_adduct")) {
    if (!col %in% names(fam_summary)) fam_summary[[col]] <- NA_real_
  }
  if (!"has_role_conflict" %in% names(fam_summary)) {
    fam_summary$has_role_conflict <- NA
  }

  fam_mass_all <- as.numeric(fam_summary$neutral_mass_consensus)
  keep_fam <- is.finite(fam_mass_all)
  fam_summary_f <- fam_summary[keep_fam, , drop = FALSE]

  family_clusters <- data.frame(
    neutral_mass_id = integer(), neutral_mass_consensus = numeric(),
    mean_score_adduct = numeric(), median_score_adduct = numeric(),
    has_role_conflict = logical()
  )
  neutral_mass_id_map <- data.frame(family_id = integer(), neutral_mass_id = integer())

  if (nrow(fam_summary_f) > 0) {
    fam_summary_f$neutral_mass_id <- assign_neutral_mass_clusters(
      fam_summary_f$neutral_mass_consensus, neutral_cluster_ppm
    )
    neutral_mass_id_map <- fam_summary_f[, c("family_id", "neutral_mass_id")]

    family_clusters <- fam_summary_f |>
      dplyr::group_by(neutral_mass_id) |>
      dplyr::summarise(
        neutral_mass_consensus = mean(neutral_mass_consensus, na.rm = TRUE),
        mean_score_adduct = mean(mean_score_adduct, na.rm = TRUE),
        median_score_adduct = mean(median_score_adduct, na.rm = TRUE),
        has_role_conflict = any(has_role_conflict, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # --- Family-side feature/evidence aggregation ---------------------------
  member_info <- fam_members |>
    dplyr::left_join(neutral_mass_id_map, by = "family_id") |>
    dplyr::filter(!is.na(neutral_mass_id))

  feature_cols <- intersect(
    c("idx", "is_c13_m0", "c13_score", "has_c13_m2_support",
      "c13_m2_score", "has_eips", "eips_elements", "eips_score"),
    names(feature_summary)
  )

  member_feature_info <- member_info |>
    dplyr::left_join(
      feature_summary |> dplyr::select(dplyr::all_of(feature_cols)),
      by = "idx"
    ) |>
    dplyr::mutate(
      is_c13_m0 = dplyr::coalesce(is_c13_m0, FALSE),
      has_c13_m2_support = dplyr::coalesce(has_c13_m2_support, FALSE),
      has_eips = dplyr::coalesce(has_eips, FALSE)
    )

  family_evidence <- member_feature_info |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      feature_idx = paste(sort(unique(idx)), collapse = ";"),
      inferred_adduct = paste(sort(unique(adduct)), collapse = ";"),
      .groups = "drop"
    )

  # --- Family candidate identity retrieval ---------------------------------
  neutral_for_matching_fam <- family_clusters |>
    dplyr::select(neutral_mass_id, neutral_mass_consensus) |>
    dplyr::distinct()

  compound_candidates_fam <- match_compound_db_by_mass(
    neutral_for_matching_fam, compound_db, ppm_tol = ppm_tol, top_n = top_n
  )

  standard_candidates_fam <- NULL
  if (use_standards && nrow(neutral_for_matching_fam) > 0) {
    standard_candidates_fam <- match_standards_by_mass(
      neutral_for_matching_fam, standards_db,
      ion_mode = ion_mode, matrix = matrix, ppm_tol = ppm_tol, quiet = quiet
    )
  }

  resolved_fam <- resolve_cross_source_identity(compound_candidates_fam, standard_candidates_fam)

  family_annotations <- resolved_fam |>
    dplyr::left_join(family_evidence, by = "neutral_mass_id") |>
    dplyr::left_join(
      family_clusters |>
        dplyr::select(
          neutral_mass_id, neutral_mass_consensus,
          mean_score_adduct, median_score_adduct, has_role_conflict
        ),
      by = "neutral_mass_id"
    ) |>
    dplyr::mutate(
      hypothesis_origin = "family",
      adduct_spatial_score = mean_score_adduct,
      family_coherence_raw = median_score_adduct,
      family_has_role_conflict = has_role_conflict
    ) |>
    dplyr::select(-mean_score_adduct, -median_score_adduct, -has_role_conflict)

  # --- Single-adduct hypotheses --------------------------------------------
  single_rows <- NULL
  if (isTRUE(include_single_adduct) && !is.null(pkm)) {
    single_rows <- build_single_adduct_candidates(
      pkm = pkm,
      feature_summary = feature_summary,
      adduct_fam = adduct_fam,
      ion_mode = ion_mode,
      matrix = matrix,
      adducts = adducts,
      compound_db = compound_db,
      standards_db = standards_db,
      ppm_tol = ppm_tol,
      neutral_cluster_ppm = neutral_cluster_ppm,
      top_n = top_n,
      quiet = quiet
    )
  }

  single_annotations <- NULL
  if (!is.null(single_rows) && nrow(single_rows) > 0) {
    single_annotations <- single_rows |>
      dplyr::mutate(feature_idx = as.character(feature_idx)) |>
      dplyr::left_join(
        family_clusters |> dplyr::select(neutral_mass_id, median_score_adduct, has_role_conflict),
        by = "neutral_mass_id"
      ) |>
      dplyr::rename(
        family_coherence_raw = median_score_adduct,
        family_has_role_conflict = has_role_conflict
      )
  }

  combined <- dplyr::bind_rows(family_annotations, single_annotations)

  if (!nrow(combined)) {
    return(empty_annotations())
  }

  # --- Neutral-mass-level isotope/EIPS evidence (family members + single
  # -adduct features together) ----------------------------------------------
  single_feature_rows <- NULL
  if (!is.null(single_rows) && nrow(single_rows) > 0) {
    single_feature_rows <- single_rows |>
      dplyr::distinct(neutral_mass_id, feature_idx) |>
      dplyr::rename(idx = feature_idx) |>
      dplyr::left_join(
        feature_summary |> dplyr::select(dplyr::all_of(feature_cols)),
        by = "idx"
      ) |>
      dplyr::mutate(
        is_c13_m0 = dplyr::coalesce(is_c13_m0, FALSE),
        has_eips = dplyr::coalesce(has_eips, FALSE)
      )
  }

  feature_evidence_cols <- c("neutral_mass_id", "idx", "is_c13_m0", "c13_score", "has_eips", "eips_score")

  all_feature_rows <- dplyr::bind_rows(
    member_feature_info |> dplyr::select(dplyr::any_of(feature_evidence_cols)),
    if (!is.null(single_feature_rows)) {
      single_feature_rows |> dplyr::select(dplyr::any_of(feature_evidence_cols))
    }
  )

  if (nrow(all_feature_rows) > 0) {
    all_feature_rows <- all_feature_rows |>
      dplyr::distinct(neutral_mass_id, idx, .keep_all = TRUE)

    evidence_by_nmid <- all_feature_rows |>
      dplyr::group_by(neutral_mass_id) |>
      dplyr::summarise(
        is_c13_m0_any = any(is_c13_m0, na.rm = TRUE),
        isotope_evidence_score = if (any(is.finite(c13_score))) max(c13_score, na.rm = TRUE) else NA_real_,
        has_eips_any = any(has_eips, na.rm = TRUE),
        eips_evidence_score = if (any(is.finite(eips_score))) max(eips_score, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
  } else {
    evidence_by_nmid <- data.frame(
      neutral_mass_id = integer(), is_c13_m0_any = logical(),
      isotope_evidence_score = numeric(), has_eips_any = logical(),
      eips_evidence_score = numeric()
    )
  }

  # --- Detected-adducts set (family members + all single-adduct
  # hypotheses sharing a neutral_mass_id) -----------------------------------
  family_adducts_nm <- member_info |> dplyr::select(neutral_mass_id, adduct)
  single_adducts_nm <- if (!is.null(single_rows) && nrow(single_rows) > 0) {
    single_rows |>
      dplyr::distinct(neutral_mass_id, inferred_adduct) |>
      dplyr::rename(adduct = inferred_adduct)
  } else {
    NULL
  }
  detected_map <- compute_detected_adducts_map(family_adducts_nm, single_adducts_nm)

  combined <- combined |>
    dplyr::left_join(evidence_by_nmid, by = "neutral_mass_id") |>
    dplyr::mutate(
      is_c13_m0 = dplyr::coalesce(is_c13_m0_any, FALSE),
      has_eips = dplyr::coalesce(has_eips_any, FALSE),
      mass_error_score = ppm_error_to_score(candidate_ppm_error, ppm_tol),
      family_coherence_score = family_coherence_raw *
        ifelse(dplyr::coalesce(family_has_role_conflict, FALSE), 0.7, 1.0)
    )

  combined$standard_adduct_recovery_score <- mapply(
    function(nmid, std_adducts) {
      if (is.na(std_adducts)) return(NA_real_)
      det <- detected_map[[as.character(nmid)]]
      if (is.null(det)) det <- character(0)
      compute_standard_adduct_recovery(std_adducts, det)
    },
    combined$neutral_mass_id, combined$standard_db_adducts
  )

  combined$priority_score <- compute_priority_score(
    combined[, c(
      "mass_error_score", "adduct_spatial_score", "isotope_evidence_score",
      "eips_evidence_score", "standard_adduct_recovery_score", "family_coherence_score"
    )],
    weights
  )

  combined <- combined |>
    dplyr::mutate(
      .broad_evidence_count =
        as.integer(is_c13_m0) + as.integer(has_eips) +
        as.integer(dplyr::coalesce(adduct_spatial_score, 0) >= adduct_min_score),
      confidence_class = dplyr::case_when(
        is.na(source) ~ NA_character_,
        source == "both" &
          dplyr::coalesce(standard_adduct_recovery_score, 0) >= recovery_threshold ~
          "cross_db_linked_adduct_recovered",
        source == "both" ~ "cross_db_linked_partial_recovery",
        source == "standards_only" ~ "standards_only",
        source == "broad_db_only" & .broad_evidence_count >= 2 ~
          "broad_db_only_multi_evidence",
        source == "broad_db_only" & .broad_evidence_count == 1 ~
          "broad_db_only_single_evidence",
        source == "broad_db_only" ~ "broad_db_only_mass_only",
        TRUE ~ NA_character_
      )
    )

  ambig <- combined |>
    dplyr::filter(!is.na(priority_score)) |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      .top1 = sort(priority_score, decreasing = TRUE)[1],
      .top2 = if (dplyr::n() >= 2) sort(priority_score, decreasing = TRUE)[2] else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ambiguous_isomeric = is.finite(.top1) & is.finite(.top2) & (.top1 - .top2) < ambiguity_gap
    ) |>
    dplyr::select(neutral_mass_id, ambiguous_isomeric)

  combined <- combined |>
    dplyr::left_join(ambig, by = "neutral_mass_id") |>
    dplyr::mutate(ambiguous_isomeric = dplyr::coalesce(ambiguous_isomeric, FALSE))

  combined |>
    dplyr::select(
      neutral_mass_id, neutral_mass_consensus, feature_idx, inferred_adduct,
      source, identity_match_type, possible_link_formula_only, hypothesis_origin,
      candidate_neutral_mass, candidate_ppm_error,
      broad_db_name, broad_db_formula, broad_db_inchikey, broad_db_inchi,
      broad_db_smiles, broad_db_id, broad_db_source,
      standard_db_name, standard_db_formula, standard_db_inchikey,
      standard_db_smiles, standard_db_compound_id, standard_db_adducts,
      mass_error_score, adduct_spatial_score, isotope_evidence_score,
      eips_evidence_score, standard_adduct_recovery_score, family_coherence_score,
      priority_score, confidence_class, ambiguous_isomeric
    ) |>
    dplyr::arrange(neutral_mass_id, dplyr::desc(priority_score))
}
