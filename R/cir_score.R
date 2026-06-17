#' CIR score
#' #' Validate C13 M+1 candidates using carbon isotope ratios
#'
#' Filters C13 M+1 candidates from `iso_morphology_candidates()` and compares
#' the observed M+1/M+0 intensity ratio against the theoretical carbon isotope
#' ratio expected for the monoisotopic m/z.
#'
#' The function is intended to add isotope-ratio evidence to morphology-based
#' C13 candidates. It does not assign definitive isotope annotations by itself.
#'
#' Filters C13_M1 candidates and validates against theoretical ratios.
#'
#' @param result Isotope detection results (`iso_morphology_candidates()` output)
#' @param pkm Peak matrix list (`mass`, `intensity`)
#' @param min_score_final Minimum `score_final` threshold in the range 0 to 1. Default: 0.6
#' @param cir_rel_tol CIR relative tolerance in the range 0 to 1. Default: 0.3
#' @param ratio_method `"sum"` (default), `"median"`, or `"mean"` for I_M1/I_M0
#' @param min_quantile Quantile mask for M0 intensities. Default 0.01
#' @param mask_strategy `"M0_only"` (recommended) or `"AND"`. Default `"M0_only"`.
#'
#' @return A data.frame with one row per C13 M+1 candidate and the following
#'   additional columns: `R_obs`, `R_theo`, `cir_rel_err`,
#'   `is_valid_c13_raw`, `cir_score`, `cir_class`, `is_chained_c13`,
#'   `is_valid_c13`, `has_C13_M2`, `idx_C13_M2`, `mz_C13_M2` and
#'   `score_C13_M2`.
#'
#' @export
cir_score <- function(result, pkm,
                      min_score_final = 0.6,
                      cir_rel_tol = 0.3,
                      ratio_method = c("sum", "mean", "median"),
                      min_quantile = 0.01,
                      mask_strategy = c("M0_only", "AND")) {

  ratio_method   <- match.arg(ratio_method)
  mask_strategy  <- match.arg(mask_strategy)

  gaussian_score <- function(err, tol) {
    err <- pmax(0, err)
    exp(-(err^2) / (2 * tol^2))
  }

  #Checks
  if (!is.data.frame(result)) {
    stop("`result` must be a data.frame returned by `iso_morphology_candidates()`.", call. = FALSE)
  }

  if (!is.list(pkm) || is.null(pkm$mass) || is.null(pkm$intensity)) {
    stop("`pkm` must be a list containing `mass` and `intensity`.", call. = FALSE)
  }

  if (!is.numeric(pkm$mass)) {
    stop("`pkm$mass` must be a numeric vector.", call. = FALSE)
  }

  if (!is.matrix(pkm$intensity) || !is.numeric(pkm$intensity)) {
    stop("`pkm$intensity` must be a numeric matrix.", call. = FALSE)
  }

  if (length(pkm$mass) != ncol(pkm$intensity)) {
    stop("Length of `pkm$mass` must match the number of columns in `pkm$intensity`.", call. = FALSE)
  }

  required_cols <- c("iso_type", "score_final", "idx_M0", "idx_cand")
  missing_cols <- setdiff(required_cols, names(result))

  if (length(missing_cols) > 0) {
    stop(
      "`result` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }


  # Filter ONLY C13_M1 with good morphology
  c13_pairs <- result |>
    dplyr::filter(iso_type == "C13_M1",
                  !is.na(score_final),
                  score_final >= min_score_final)

  if (nrow(c13_pairs) == 0) {
    message("No C13_M1 candidates with score_final >= ", min_score_final)
    return(c13_pairs)
  }

  # Indices
  idx_M0 <- c13_pairs$idx_M0
  idx_M1 <- c13_pairs$idx_cand

  # Intensities
  I_M0 <- pkm$intensity[, idx_M0, drop = FALSE]
  I_M1 <- pkm$intensity[, idx_M1, drop = FALSE]

  # Quantile thresholds (per feature)
  min_I_M0 <- apply(I_M0, 2, function(x) stats::quantile(x, probs = min_quantile, na.rm = TRUE))
  min_I_M1 <- apply(I_M1, 2, function(x) stats::quantile(x, probs = min_quantile, na.rm = TRUE))

  # Mask pixels (recommended: M0_only for C13)
  if (mask_strategy == "M0_only") {
    mask_valid <- sweep(I_M0, 2, min_I_M0, FUN = ">")
  } else {
    mask_valid <- sweep(I_M0, 2, min_I_M0, FUN = ">") & sweep(I_M1, 2, min_I_M1, FUN = ">")
  }

  # Apply the mask
  I_M0_f <- I_M0
  I_M1_f <- I_M1
  I_M0_f[!mask_valid] <- NA
  I_M1_f[!mask_valid] <- NA

  # Observed ratios
  R_obs_pixel <- I_M1_f / pmax(I_M0_f, 1e-6)

  R_obs <- if (ratio_method == "mean") {
    colMeans(R_obs_pixel, na.rm = TRUE)
  } else if (ratio_method == "median") {
    apply(R_obs_pixel, 2, stats::median, na.rm = TRUE)
  } else {
    # IMPORTANT: use FILTERED intensities (I_M*_f), not raw
    colSums(I_M1_f, na.rm = TRUE) / pmax(colSums(I_M0_f, na.rm = TRUE), 1e-6)
  }

  # Theoretical CIR (M+1 only)
  mz_mono <- pkm$mass[idx_M0]
  R_theo  <- cir_ratio(mz_mono)


  # CIR validation
  cir_rel_err <- abs(R_obs - R_theo) / pmax(R_theo, 1e-6)

  is_valid_c13_raw <- is.finite(cir_rel_err) & cir_rel_err < cir_rel_tol

  cir_score <- gaussian_score(cir_rel_err, cir_rel_tol)

  cir_class <- dplyr::case_when(
    is_valid_c13_raw ~ "high_agreement",
    cir_rel_err < 1.5 * cir_rel_tol ~ "moderate_agreement",
    TRUE ~ "low_agreement"
  )

  # Build output table
  out <- data.frame(
    idx_M0           = idx_M0,
    mz_M0            = mz_mono,
    idx_M1           = idx_M1,
    mz_M1            = pkm$mass[idx_M1],
    R_obs            = as.numeric(R_obs),
    R_theo           = as.numeric(R_theo),
    cir_rel_err      = as.numeric(cir_rel_err),
    is_valid_c13_raw = as.logical(is_valid_c13_raw),
    cir_score        = as.numeric(cir_score),
    cir_class        = as.character(cir_class),
    stringsAsFactors = FALSE
  )

  # Detect chained isotope artefacts:
  # if a peak is already a valid M+1,
  # it should not be reused as a new M0
  valid_m1 <- out$idx_M1[out$is_valid_c13_raw %in% TRUE]

  out$is_chained_c13 <- out$idx_M0 %in% valid_m1

  # Final validity
  out$is_valid_c13 <- out$is_valid_c13_raw & !out$is_chained_c13
  # Add C13 M+2 support for CIR-valid M0 features
  c13_m2 <- result |>
    dplyr::filter(
      iso_type == "C13_M2",
      !is.na(score_final),
      score_final >= min_score_final
    ) |>
    dplyr::group_by(idx_M0) |>
    dplyr::slice_max(order_by = score_final, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx_M0,
      has_C13_M2 = TRUE,
      idx_C13_M2 = idx_cand,
      mz_C13_M2 = pkm$mass[idx_cand],
      score_C13_M2 = score_final
    )

  out <- out |>
    dplyr::left_join(c13_m2, by = "idx_M0") |>
    dplyr::mutate(
      has_C13_M2 = dplyr::coalesce(has_C13_M2, FALSE),
      idx_C13_M2 = dplyr::if_else(
        is_valid_c13 & has_C13_M2,
        as.integer(idx_C13_M2),
        NA_integer_
      ),
      mz_C13_M2 = dplyr::if_else(
        is_valid_c13 & has_C13_M2,
        as.numeric(mz_C13_M2),
        NA_real_
      ),
      score_C13_M2 = dplyr::if_else(
        is_valid_c13 & has_C13_M2,
        as.numeric(score_C13_M2),
        NA_real_
      ),
      has_C13_M2 = is_valid_c13 & has_C13_M2
    )
  out
}

# Example:
# %>% %>% cir_results <- cir_score(res_morph, pkm, min_score_final = 0.6, ratio_method = "mean")
