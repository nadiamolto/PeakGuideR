#' Run the PeakGuideR annotation workflow
#'
#' @description
#' Runs the main PeakGuideR workflow from either a PeakGuideR peak matrix object
#' or a Cardinal MSI object.
#'
#' The workflow includes isotope morphology detection, carbon isotope-ratio
#' validation, elemental isotope-pattern support, adduct candidate detection,
#' adduct-family grouping, relation-table construction, feature-level
#' summarization and neutral-mass candidate matching.
#'
#' @param pkm rMSI2 peak matrix object, or a supported Cardinal
#'   `MSImagingExperiment` object. Cardinal objects are converted internally
#'   using `cardinal_to_peakmatrix()`.
#' @param ion_mode Ion mode, either `"pos"` or `"neg"`.
#' @param matrix Matrix name. Use `"HCCA"` to enable HCCA-specific standard
#'   adduct support in neutral-mass candidate matching.
#' @param morph_prefer_mode Mass-deviation preference mode used by
#'   `iso_morphology_candidates()`, either `"ppm"` or `"dp"`.
#' @param morph_method Spatial similarity method used by
#'   `iso_morphology_candidates()`.
#' @param morph_transform Intensity transformation used by
#'   `iso_morphology_candidates()`.
#' @param morph_tile_blend Method used to combine tile-level morphology scores.
#' @param adducts Optional adduct definition table. If `NULL`, PeakGuideR uses
#'   `default_adducts(ion_mode)`. Users can inspect and modify the default
#'   adduct table with `default_adducts()`.
#' @param compound_db Optional compound mass database.
#' @param standards_db Optional standard adduct library.
#' @param iso_min_score Minimum isotope morphology score used for CIR/EIPS input.
#' @param cir_rel_tol Relative tolerance for carbon isotope-ratio validation.
#' @param eips_rel_tol Relative tolerance for EIPS validation.
#' @param ratio_method Ratio aggregation method.
#' @param adduct_tol_ppm PPM tolerance for adduct candidate detection.
#' @param adduct_neutral_tol_ppm PPM tolerance for adduct neutral-mass consistency.
#' @param adduct_method Spatial similarity method used for adduct detection.
#' @param adduct_transform Intensity transformation used for adduct detection.
#' @param adduct_min_quantile Minimum quantile used for adduct spatial vectors.
#' @param adduct_clip_negatives Logical. If `TRUE`, negative transformed values
#'   are clipped in adduct detection.
#' @param adduct_min_score Minimum spatial score for adduct candidates/families.
#' @param neutral_cluster_ppm PPM tolerance used to cluster inferred neutral masses.
#' @param candidate_ppm_tol PPM tolerance for compound candidate matching.
#' @param top_n Maximum number of compound candidates per neutral mass.
#' @param quiet Logical. If `FALSE`, prints progress messages.
#'
#' @return A list with all main PeakGuideR workflow outputs.
#' @export
#'
run_peakguider_workflow <- function(
    pkm,
    ion_mode = c("pos", "neg"),
    matrix = NULL, # Standard-adduct support is currently matrix-specific. At present, the included standard library is available only for HCCA+DEA matrix.
    adducts = NULL,
    compound_db = NULL,
    standards_db = NULL,
    morph_prefer_mode = c("ppm", "dp"),
    morph_method = c("pearson", "cosine", "spearman"),
    morph_transform = c("none", "log1p", "zscore"),
    morph_tile_blend = c("median", "p25", "pass_rate"),
    iso_min_score = 0.6,
    cir_rel_tol = 0.3,
    eips_rel_tol = 0.3,
    ratio_method = c("sum", "mean", "median"),
    adduct_tol_ppm = 5,
    adduct_neutral_tol_ppm = 5,
    adduct_method = c("pearson", "cosine", "spearman"),
    adduct_transform = c("none", "log1p", "zscore"),
    adduct_min_quantile = 0.01,
    adduct_clip_negatives = TRUE,
    adduct_min_score = 0.5,
    neutral_cluster_ppm = 5,
    candidate_ppm_tol = 5,
    top_n = 10L,
    quiet = FALSE
) {

  ion_mode <- match.arg(ion_mode)
  morph_prefer_mode <- match.arg(morph_prefer_mode)
  morph_method <- match.arg(morph_method)
  morph_transform <- match.arg(morph_transform)
  morph_tile_blend <- match.arg(morph_tile_blend)
  ratio_method <- match.arg(ratio_method)
  adduct_method <- match.arg(adduct_method)
  adduct_transform <- match.arg(adduct_transform)


  input_type <- "peak_matrix"

  if (is.null(adducts)) {
    adducts <- default_adducts(ion_mode)

    if (!isTRUE(quiet)) {
      message(
        "No adduct table supplied; using PeakGuideR default adduct definitions for ",
        ion_mode,
        " mode."
      )
    }}

  if (missing(pkm)) {
    stop(
      "Argument 'pkm' is required. Please provide a peak matrix or a supported Cardinal object.",
      call. = FALSE
    )
  }

  if (is_cardinal_object(pkm)) {
    input_type <- "cardinal_object"

    if (!isTRUE(quiet)) {
      message("Converting Cardinal object to PeakGuideR peak matrix...")
    }

    pkm <- cardinal_to_peakmatrix(pkm)
  }

  if (!is.list(pkm) || is.null(pkm$mass)) {
    stop(
      "`pkm` must be a peak matrix object containing at least `mass`, `intensity` and `pos`, or a supported Cardinal object.",
      call. = FALSE
    )
  }

  if (!isTRUE(quiet)) {
    message("1/8 Detecting isotope morphology candidates...")
  }

  morph_results <- iso_morphology_candidates(
    pkm = pkm,
    prefer_mode = morph_prefer_mode,
    method = morph_method,
    transform = morph_transform,
    tile_blend = morph_tile_blend
  )

  if (!isTRUE(quiet)) {
    message("2/8 Running CIR validation...")
  }

  cir_results <- cir_score(
    result = morph_results,
    pkm = pkm,
    min_score_final = iso_min_score,
    cir_rel_tol = cir_rel_tol,
    ratio_method = ratio_method
  )

  if (!isTRUE(quiet)) {
    message("3/8 Running EIPS validation...")
  }


  # EIPS reference tables are internal package data; they are not exposed in the main workflow interface.

  eips_n_ref <- get("eips_n_table", envir = asNamespace("PeakGuideR"))
  eips_ref <- get("eips_table", envir = asNamespace("PeakGuideR"))


  eips_results <- eips_score(
    result = morph_results,
    pkm = pkm,
    eips_n_table = eips_n_ref,
    eips_table = eips_ref,
    ion_mode = ion_mode,
    adducts = adducts,
    min_score_final = iso_min_score,
    eips_rel_tol = eips_rel_tol,
    ratio_method = ratio_method,
    cir_df = cir_results,
    morph_df = morph_results
  )

  if (!isTRUE(quiet)) {
    message("4/8 Detecting adduct candidates...")
  }

  adduct_edges <- adduct_candidates(
    pkm = pkm,
    ion_mode = ion_mode,
    adducts = adducts,
    tol_ppm = adduct_tol_ppm,
    neutral_tol_ppm = adduct_neutral_tol_ppm,
    method = adduct_method,
    transform = adduct_transform,
    min_quantile = adduct_min_quantile,
    clip_negatives = adduct_clip_negatives,
    min_score_spatial = adduct_min_score
  )

  if (!isTRUE(quiet)) {
    message("5/8 Building adduct families...")
  }

  adduct_fam <- adduct_families(
    adduct_edges = adduct_edges,
    min_score_adduct = adduct_min_score,
    neutral_cluster_ppm = neutral_cluster_ppm,
    min_family_size = 2L,
    use_only_valid = TRUE
  )

  if (!isTRUE(quiet)) {
    message("6/8 Building relation table...")
  }

  relation_table <- build_relation_table(
    cir_results = cir_results,
    eips_results = eips_results,
    adduct_fam = adduct_fam,
    only_valid = TRUE
  )

  if (!isTRUE(quiet)) {
    message("7/8 Building feature summary...")
  }

  feature_summary <- build_feature_summary(
    relation_table = relation_table,
    adduct_fam = adduct_fam,
    pkm = pkm
  )

  if (!isTRUE(quiet)) {
    message("8/8 Building neutral-mass candidates...")
  }

  neutral_mass_candidates <- build_neutral_mass_candidates(
    adduct_fam = adduct_fam,
    feature_summary = feature_summary,
    compound_db = compound_db,
    standards_db = standards_db,
    ion_mode = ion_mode,
    matrix = matrix,
    ppm_tol = candidate_ppm_tol,
    neutral_cluster_ppm = neutral_cluster_ppm,
    top_n = top_n,
    quiet = quiet
  )

  out <- list(
    morph_results = morph_results,
    cir_results = cir_results,
    eips_results = eips_results,
    adduct_edges = adduct_edges,
    adduct_families = adduct_fam,
    relation_table = relation_table,
    feature_summary = feature_summary,
    neutral_mass_candidates = neutral_mass_candidates,
    pkm = pkm,
    parameters = list(
      input_type = input_type,
      ion_mode = ion_mode,
      matrix = matrix,
      morph_prefer_mode = morph_prefer_mode,
      morph_method = morph_method,
      morph_transform = morph_transform,
      morph_tile_blend = morph_tile_blend,
      iso_min_score = iso_min_score,
      cir_rel_tol = cir_rel_tol,
      eips_rel_tol = eips_rel_tol,
      ratio_method = ratio_method,
      adduct_tol_ppm = adduct_tol_ppm,
      adduct_neutral_tol_ppm = adduct_neutral_tol_ppm,
      adduct_method = adduct_method,
      adduct_transform = adduct_transform,
      adduct_min_quantile = adduct_min_quantile,
      adduct_clip_negatives = adduct_clip_negatives,
      adduct_min_score = adduct_min_score,
      neutral_cluster_ppm = neutral_cluster_ppm,
      candidate_ppm_tol = candidate_ppm_tol,
      top_n = top_n
    )
  )

  class(out) <- c("peakguider_workflow", class(out))

  out
}


#' Print PeakGuideR workflow result
#'
#' @param x A `peakguider_workflow` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.peakguider_workflow <- function(x, ...) {
  cat("PeakGuideR workflow result\n")
  cat("---------------------------\n")
  cat("Input type:                 ", x$parameters$input_type, "\n")
  cat("Ion mode:                   ", x$parameters$ion_mode, "\n")
  cat("Matrix:                     ", ifelse(is.null(x$parameters$matrix), "none", x$parameters$matrix), "\n")
  cat("Morphology candidates:      ", nrow(x$morph_results), "\n")
  cat("CIR results:                ", nrow(x$cir_results), "\n")
  cat("EIPS results:               ", nrow(x$eips_results), "\n")
  cat("Adduct candidate edges:     ", nrow(x$adduct_edges), "\n")

  if (is.list(x$adduct_families) &&
      "family_summary" %in% names(x$adduct_families)) {
    cat("Adduct families:            ", nrow(x$adduct_families$family_summary), "\n")
  } else {
    cat("Adduct families:             NA\n")
  }

  cat("Relation table:             ", nrow(x$relation_table), "\n")
  cat("Feature summary:            ", nrow(x$feature_summary), "\n")
  cat("Neutral-mass candidates:    ", nrow(x$neutral_mass_candidates), "\n")

  invisible(x)
}
