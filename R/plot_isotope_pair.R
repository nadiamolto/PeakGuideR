#' Plot an isotope pair (and, when available, its M+2 satellite)
#'
#' @description
#' Plots the ion image of a candidate monoisotopic feature (M0) side by side
#' with its candidate isotope satellite (M+k), reusing `plot_ion_image()` for
#' each panel. If a `cir_score()` row is supplied and it has a supported M+2
#' satellite (`has_C13_M2 == TRUE`), a third panel for that M+2 feature is
#' added automatically.
#'
#' @param pkm Peak matrix list, see `plot_ion_image()`.
#' @param idx_M0 Integer. Column index of the monoisotopic (M0) feature.
#' @param idx_cand Integer. Column index of the candidate isotope satellite
#'   feature (M+k).
#' @param morph_row Optional single-row data.frame - one row of the output of
#'   `iso_morphology_candidates()` - used (if `cir_row` and `eips_row` are not
#'   supplied) to annotate the plot subtitle with `score_final`.
#' @param cir_row Optional single-row data.frame - one row of the output of
#'   `cir_score()`. When supplied, its `R_obs`/`R_theo`/`cir_class` are used
#'   for the plot subtitle (taking priority over `eips_row`/`morph_row`), and
#'   an M+2 panel is added automatically when `has_C13_M2` is `TRUE`.
#' @param eips_row Optional single-row data.frame - one row of the output of
#'   `eips_score()`. Used for the plot subtitle when `cir_row` is not
#'   supplied.
#' @param ... Additional arguments passed through to `plot_ion_image()` (e.g.
#'   `palette`, `clip_quantile`, `flip_y`, `show_legend`).
#'
#' @return A `patchwork` object combining the M0/M+k (and, when available,
#'   M+2) ion image panels.
#'
#' @examples
#' \dontrun{
#' data(example_pkm, package = "PeakGuideR")
#'
#' morph_results <- iso_morphology_candidates(example_pkm, prefer_mode = "ppm")
#' cir_results <- cir_score(morph_results, example_pkm)
#'
#' one_row <- cir_results[1, ]
#' plot_isotope_pair(
#'   example_pkm,
#'   idx_M0 = one_row$idx_M0,
#'   idx_cand = one_row$idx_M1,
#'   cir_row = one_row
#' )
#' }
#'
#' @export
plot_isotope_pair <- function(
    pkm,
    idx_M0,
    idx_cand,
    morph_row = NULL,
    cir_row = NULL,
    eips_row = NULL,
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required to use `plot_isotope_pair()`. ",
      "Please install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required to use `plot_isotope_pair()`. ",
      "Please install it with install.packages('patchwork').",
      call. = FALSE
    )
  }

  if (!is.null(morph_row)) {
    if (!is.data.frame(morph_row) || nrow(morph_row) != 1) {
      stop("`morph_row` must be a single-row data.frame from `iso_morphology_candidates()`.", call. = FALSE)
    }
    required_cols <- c("idx_M0", "mz_M0", "idx_cand", "mz_cand", "score_final", "iso_type", "element")
    missing_cols <- setdiff(required_cols, names(morph_row))
    if (length(missing_cols) > 0) {
      stop(
        "`morph_row` is missing required columns: ", paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(cir_row)) {
    if (!is.data.frame(cir_row) || nrow(cir_row) != 1) {
      stop("`cir_row` must be a single-row data.frame from `cir_score()`.", call. = FALSE)
    }
    required_cols <- c(
      "idx_M0", "mz_M0", "idx_M1", "mz_M1", "R_obs", "R_theo", "cir_score", "cir_class",
      "is_valid_c13", "has_C13_M2", "idx_C13_M2", "mz_C13_M2", "score_C13_M2"
    )
    missing_cols <- setdiff(required_cols, names(cir_row))
    if (length(missing_cols) > 0) {
      stop(
        "`cir_row` is missing required columns: ", paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(eips_row)) {
    if (!is.data.frame(eips_row) || nrow(eips_row) != 1) {
      stop("`eips_row` must be a single-row data.frame from `eips_score()`.", call. = FALSE)
    }
    required_cols <- c(
      "idx_M0", "mz_M0", "idx_cand", "mz_iso", "element",
      "R_obs", "n_hat", "R_theo_hat", "score_eips", "is_valid_eips"
    )
    missing_cols <- setdiff(required_cols, names(eips_row))
    if (length(missing_cols) > 0) {
      stop(
        "`eips_row` is missing required columns: ", paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  p_m0 <- plot_ion_image(
    pkm, idx = idx_M0,
    title = sprintf("M0 - m/z %.4f", pkm$mass[idx_M0]),
    ...
  )
  p_cand <- plot_ion_image(
    pkm, idx = idx_cand,
    title = sprintf("M+k - m/z %.4f", pkm$mass[idx_cand]),
    ...
  )

  panels <- p_m0 | p_cand

  if (!is.null(cir_row) && isTRUE(cir_row$has_C13_M2)) {
    p_m2 <- plot_ion_image(
      pkm, idx = cir_row$idx_C13_M2,
      title = sprintf("M+2 - m/z %.4f", cir_row$mz_C13_M2),
      ...
    )
    panels <- panels | p_m2
  }

  subtitle <- if (!is.null(cir_row)) {
    sprintf(
      "CIR: R_obs=%.3f vs R_theo=%.3f (%s)",
      cir_row$R_obs, cir_row$R_theo, cir_row$cir_class
    )
  } else if (!is.null(eips_row)) {
    sprintf("EIPS (%s): score=%.3f", eips_row$element, eips_row$score_eips)
  } else if (!is.null(morph_row)) {
    sprintf("morphology score_final=%.3f", morph_row$score_final)
  } else {
    NULL
  }

  panels + patchwork::plot_annotation(subtitle = subtitle)
}
