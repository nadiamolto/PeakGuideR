#' Build neutral-mass candidate table
#'
#' @description
#' Builds a neutral-mass candidate table from adduct families.
#'
#' This function is a per-`neutral_mass_id` summary view derived from
#' [build_candidate_annotations()] (family-derived candidates only; see
#' `include_single_adduct` there for the full single-adduct-aware table
#' returned as `candidate_annotations` by [run_peakguider_workflow()]).
#' Compound and standard-adduct identity are resolved with
#' [resolve_cross_source_identity()], so a shared molecular formula alone
#' never fuses two different compounds into the same
#' `has_standard_compound_match` (see the `dontrun` example in
#' [resolve_cross_source_identity()] for details).
#'
#' The output contains one row per inferred neutral mass and compound/standard
#' candidate. Candidate compounds are putative mass matches, not definitive
#' identifications.
#'
#' @param adduct_fam Output from `adduct_families()`.
#' @param feature_summary Output from `build_feature_summary()`.
#' @param compound_db Compound mass database. If `NULL`, the included
#'   non-commercial compound mass database is loaded.
#' @param standards_db Standard adduct library. If `NULL` and `matrix = "HCCA"`,
#'   the included non-commercial standard adduct library is loaded.
#' @param ion_mode `"pos"` or `"neg"`.
#' @param matrix Matrix name. Standard-adduct support is currently applied only
#'   when `matrix = "HCCA"`. Use `NULL` to skip standard-adduct support.
#' @param ppm_tol PPM tolerance for compound mass matching.
#' @param neutral_cluster_ppm PPM tolerance used to group similar neutral masses
#'   inferred from different adduct families.
#' @param top_n Maximum number of compound candidates kept per neutral mass.
#' Use `NULL` to retain all candidates within `ppm_tol`.
#' @param candidate_annotations Optional pre-computed output from
#'   [build_candidate_annotations()]. If `NULL` (the default), the function
#'   computes its own family-only candidate table internally, exactly as
#'   before. If supplied (for example the same table already computed by
#'   [run_peakguider_workflow()], possibly including single-adduct rows),
#'   it is reused instead of recomputing the compound/standard matching and
#'   identity resolution from scratch; rows other than
#'   `hypothesis_origin == "family"` are filtered out first, so the returned
#'   table is unaffected by whether single-adduct hypotheses were included
#'   upstream.
#' @param quiet Logical. If `FALSE`, database loading functions may print notices.
#'
#' @return A data.frame with one row per neutral mass and compound/standard
#'   candidate.
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
#' neutral_mass_candidates <- build_neutral_mass_candidates(
#'   adduct_fam = adduct_fam,
#'   feature_summary = feature_summary,
#'   ion_mode = "pos",
#'   matrix = "HCCA"
#' )
#' head(neutral_mass_candidates)
#' }
#' @export
build_neutral_mass_candidates <- function(
    adduct_fam,
    feature_summary,
    compound_db = NULL,
    standards_db = NULL,
    ion_mode = c("pos", "neg"),
    matrix = NULL,
    ppm_tol = 5,
    neutral_cluster_ppm = 5,
    top_n = 10L,
    candidate_annotations = NULL,
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
  missing_family_cols <- setdiff(
    required_family_cols,
    names(adduct_fam$family_summary)
  )

  if (length(missing_family_cols) > 0) {
    stop(
      "adduct_fam$family_summary is missing required columns: ",
      paste(missing_family_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_member_cols <- c("family_id", "idx", "mz", "adduct")
  missing_member_cols <- setdiff(
    required_member_cols,
    names(adduct_fam$family_members)
  )

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
      stop(
        "`top_n` must be NULL or a positive integer.",
        call. = FALSE
      )
    }
  }

  ppm_tol <- as.numeric(ppm_tol)
  if (!is.finite(ppm_tol) || ppm_tol <= 0) {
    ppm_tol <- 5
  }

  neutral_cluster_ppm <- as.numeric(neutral_cluster_ppm)
  if (!is.finite(neutral_cluster_ppm) || neutral_cluster_ppm <= 0) {
    neutral_cluster_ppm <- ppm_tol
  }

  if (is.null(candidate_annotations)) {
    if (is.null(compound_db)) {
      compound_db <- load_compound_mass_database(quiet = quiet)
    }

    stopifnot(is.data.frame(compound_db))

    if (!"MonoisotopicMass" %in% names(compound_db)) {
      stop(
        "compound_db must contain a `MonoisotopicMass` column.",
        call. = FALSE
      )
    }
  } else {
    stopifnot(is.data.frame(candidate_annotations))

    if (!"hypothesis_origin" %in% names(candidate_annotations)) {
      stop(
        "`candidate_annotations` must contain a `hypothesis_origin` column.",
        call. = FALSE
      )
    }
  }

  fam_summary <- adduct_fam$family_summary
  fam_members <- adduct_fam$family_members

  fam_summary$.neutral_mass <- as.numeric(fam_summary$neutral_mass_consensus)
  keep_fam <- is.finite(fam_summary$.neutral_mass)
  fam_summary_f <- fam_summary[keep_fam, , drop = FALSE]

  if (nrow(fam_summary_f) == 0) {
    return(dplyr::tibble())
  }

  fam_summary_f$neutral_mass_id <- assign_neutral_mass_clusters(
    fam_summary_f$.neutral_mass, neutral_cluster_ppm
  )

  member_info <- fam_members |>
    dplyr::left_join(
      fam_summary_f |> dplyr::select(family_id, neutral_mass_id),
      by = "family_id"
    ) |>
    dplyr::filter(!is.na(neutral_mass_id))

  feature_cols <- intersect(
    c(
      "idx",
      "is_c13_m0",
      "c13_score",
      "has_c13_m2_support",
      "c13_m2_score",
      "has_eips",
      "eips_elements",
      "eips_score"
    ),
    names(feature_summary)
  )

  member_feature_info <- member_info |>
    dplyr::left_join(
      feature_summary |>
        dplyr::select(dplyr::all_of(feature_cols)),
      by = "idx"
    ) |>
    dplyr::mutate(
      is_c13_m0 = dplyr::coalesce(is_c13_m0, FALSE),
      has_c13_m2_support = dplyr::coalesce(has_c13_m2_support, FALSE),
      has_eips = dplyr::coalesce(has_eips, FALSE)
    )

  neutral_summary <- fam_summary_f |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      neutral_mass_consensus = mean(.neutral_mass, na.rm = TRUE),
      n_adduct_families = dplyr::n_distinct(family_id),
      adduct_family_ids = paste(sort(unique(family_id)), collapse = ";"),
      .groups = "drop"
    )

  adduct_support_summary <- member_feature_info |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      n_features = dplyr::n_distinct(idx),
      feature_idx = paste(sort(unique(idx)), collapse = ";"),
      feature_mz = paste(round(sort(unique(mz)), 6), collapse = ";"),
      n_adducts_inferred = dplyr::n_distinct(adduct),
      inferred_adducts = paste(sort(unique(adduct)), collapse = ";"),

      has_CIR_support = any(is_c13_m0, na.rm = TRUE),
      CIR_feature_idx = paste(
        sort(unique(idx[is_c13_m0 %in% TRUE])),
        collapse = ";"
      ),
      CIR_score = if (any(is.finite(c13_score))) {
        max(c13_score, na.rm = TRUE)
      } else {
        NA_real_
      },

      has_C13_M2_support = any(has_c13_m2_support, na.rm = TRUE),
      C13_M2_score = if (any(is.finite(c13_m2_score))) {
        max(c13_m2_score, na.rm = TRUE)
      } else {
        NA_real_
      },

      has_EIPS_support = any(has_eips, na.rm = TRUE),
      EIPS_elements = paste(
        sort(unique(unlist(strsplit(
          paste(eips_elements[has_eips %in% TRUE], collapse = ";"),
          ";",
          fixed = TRUE
        )))),
        collapse = ";"
      ),
      EIPS_score = if (any(is.finite(eips_score))) {
        max(eips_score, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      CIR_feature_idx = dplyr::if_else(
        CIR_feature_idx == "",
        NA_character_,
        CIR_feature_idx
      ),
      EIPS_elements = dplyr::if_else(
        EIPS_elements == "",
        NA_character_,
        EIPS_elements
      )
    )

  neutral_table <- neutral_summary |>
    dplyr::left_join(adduct_support_summary, by = "neutral_mass_id")

  # --- Candidate identity, resolved with the same identifier hierarchy
  # used across PeakGuideR (see resolve_cross_source_identity()). A
  # pre-computed table can be supplied to avoid repeating this matching when
  # the caller (e.g. run_peakguider_workflow()) already computed it. --------
  if (is.null(candidate_annotations)) {
    candidate_annotations <- build_candidate_annotations(
      adduct_fam = adduct_fam,
      feature_summary = feature_summary,
      pkm = NULL,
      ion_mode = ion_mode,
      matrix = matrix,
      compound_db = compound_db,
      standards_db = standards_db,
      ppm_tol = ppm_tol,
      neutral_cluster_ppm = neutral_cluster_ppm,
      top_n = top_n,
      include_single_adduct = FALSE,
      quiet = quiet
    )
  } else {
    candidate_annotations <- candidate_annotations |>
      dplyr::filter(hypothesis_origin == "family")
  }

  matched_inferred_adducts <- function(inferred_adduct, standard_db_adducts) {
    inferred_set <- unique(trimws(strsplit(inferred_adduct, ";", fixed = TRUE)[[1]]))
    std_set <- if (is.na(standard_db_adducts)) {
      character(0)
    } else {
      unique(trimws(strsplit(standard_db_adducts, ";", fixed = TRUE)[[1]]))
    }
    intersect(inferred_set, std_set)
  }

  candidate_view <- candidate_annotations |>
    dplyr::mutate(
      candidate_source = dplyr::if_else(
        source == "standards_only", "standards_db", broad_db_source
      ),
      candidate_db_id = dplyr::coalesce(broad_db_id, standard_db_compound_id),
      candidate_name = dplyr::coalesce(broad_db_name, standard_db_name),
      candidate_formula = dplyr::coalesce(broad_db_formula, standard_db_formula),
      candidate_inchi = broad_db_inchi,
      candidate_inchikey = dplyr::coalesce(broad_db_inchikey, standard_db_inchikey),
      candidate_smiles = dplyr::coalesce(broad_db_smiles, standard_db_smiles),
      candidate_kegg = NA_character_,
      candidate_mass = candidate_neutral_mass,
      has_standard_compound_match = source == "both",
      has_standard_adduct_match = dplyr::coalesce(standard_adduct_recovery_score, 0) > 0,
      standard_matched_name = standard_db_name,
      standard_matched_adducts = standard_db_adducts
    )

  if (nrow(candidate_view) > 0) {
    matched_list <- mapply(
      matched_inferred_adducts,
      candidate_view$inferred_adduct, candidate_view$standard_db_adducts,
      SIMPLIFY = FALSE
    )
    candidate_view$standard_matched_inferred_adducts <- vapply(
      matched_list,
      function(x) if (length(x)) paste(sort(x), collapse = ";") else NA_character_,
      character(1)
    )
    candidate_view$n_standard_matched_adducts <- lengths(matched_list)
  } else {
    candidate_view$standard_matched_inferred_adducts <- character(0)
    candidate_view$n_standard_matched_adducts <- integer(0)
  }

  candidate_view <- candidate_view |>
    dplyr::select(
      neutral_mass_id,
      candidate_source, candidate_db_id, candidate_name, candidate_formula,
      candidate_neutral_mass, candidate_inchi, candidate_inchikey,
      candidate_smiles, candidate_kegg, candidate_mass, candidate_ppm_error,
      has_standard_compound_match, has_standard_adduct_match,
      standard_matched_name, standard_matched_adducts,
      standard_matched_inferred_adducts, n_standard_matched_adducts
    )

  out <- neutral_table |>
    dplyr::left_join(candidate_view, by = "neutral_mass_id") |>
    dplyr::mutate(
      matrix = if (is.null(matrix)) NA_character_ else as.character(matrix)
    )

  out |>
    dplyr::arrange(
      neutral_mass_id,
      candidate_ppm_error
    )
}
