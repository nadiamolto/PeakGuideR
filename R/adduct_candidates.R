#' Detects adduct pairwise candidates
#'
#' @description
#' Detects pairwise adduct candidates by combining:
#' \itemize{
#'   \item expected m/z differences between adduct species
#'   \item spatial/intensity similarity between features
#'   \item agreement between inferred neutral masses.
#' }
#'
#' The function returns candidate adduct edges that can later be grouped into
#' larger ion families or pseudo-compounds.
#'
#' @param pkm A list with:
#'   \itemize{
#'     \item mass: numeric vector of m/z values.
#'     \item intensity: numeric matrix pixels by features.
#'   }
#' @param ion_mode Ionization mode. Either "pos" or "neg".
#' @param adducts Optional data-frame of adduct definitions. If NULL,
#'   default_adducts() is used (adducts usually found in MALDI-MSI).
#' @param feature_idx Optional integer vector of feature indices to evaluate.
#'   If NULL, all features are used.
#' @param tol_ppm Numeric. PPM tolerance used to search candidate peak pairs.
#'   Default 5.
#' @param neutral_tol_ppm Numeric. PPM tolerance used to score agreement
#'   between inferred neutral masses. Default `5`. This value was chosen for
#'   high-resolution Orbitrap data and should be adjusted for other instruments
#'   or preprocessing settings.
#' @param method Similarity metric between feature intensity profiles:
#'   "pearson", "cosine" or "spearman". Default "pearson".
#' @param transform Intensity transformation before similarity calculation:
#'   "none", "log1p" or "zscore". Default "none".
#' @param min_quantile Numeric in the range 0 to 1. Feature-wise low-intensity
#'   quantile filter applied after retaining co-detected pixels. Pixels are kept
#'   only when both features are above their respective threshold. Default `0.01`.
#' @param clip_negatives Logical. If TRUE, negative intensities are truncated
#'   to zero before transformation. Default TRUE.
#' @param min_score_spatial Numeric in the range 0 to 1. Minimum spatial score required
#'   to keep an adduct edge. Default 0.5.
#'
#' @details
#' The function assumes that the mass column in the adduct table stores the
#' net m/z shift relative to the neutral mass:
#'
#' \deqn{mz = M + shift}
#'
#' so the neutral mass is recovered as:
#'
#' \deqn{M = mz - shift}
#'
#' For each pair of adducts within the selected ionization mode, the expected
#' m/z difference is:
#'
#' \deqn{\Delta m/z = shift_b - shift_a}
#'
#' Candidate feature pairs matching that difference are then evaluated using
#' spatial/intensity similarity and neutral-mass consistency.
#'
#' The reported `score_adduct` corresponds to the spatial similarity score after
#' filtering candidate pairs by the expected adduct m/z difference. The
#' neutral-mass consistency score is returned as `score_mass` for inspection and
#' downstream filtering.
#'
#' @return A data.frame with one row per candidate adduct relation:
#' \itemize{
#'   \item idx_i, mz_i: first feature index and m/z
#'   \item idx_j, mz_j: second feature index and m/z
#'   \item adduct_i, adduct_j: adduct hypotheses
#'   \item delta_theo, delta_obs, delta_err_ppm
#'   \item neutral_mass_i, neutral_mass_j, neutral_mass_mean
#'   \item neutral_err_ppm
#'   \item score_spatial, score_mass, score_adduct
#'   \item is_valid_adduct
#' }
#'
#' @examples
#' \dontrun{
#' adduct_res <- adduct_candidates(
#'   pkm,
#'   ion_mode = "pos"
#' )
#'
#' extra_adducts <- rbind(
#'   default_adducts(),
#'   data.frame(
#'     name = "[M-H]-",
#'     mode = "neg",
#'     mass = -1.007276,
#'     stringsAsFactors = FALSE
#'   )
#' )
#'
#' adduct_res_neg <- adduct_candidates(
#'   pkm,
#'   ion_mode = "neg",
#'   adducts = extra_adducts
#' )
#' }
#' @export
adduct_candidates <- function(
    pkm,
    ion_mode = c("pos", "neg"),
    adducts = NULL,
    feature_idx = NULL,
    tol_ppm = 5,
    neutral_tol_ppm = 5,
    method = c("pearson", "cosine", "spearman"),
    transform = c("none", "log1p", "zscore"),
    min_quantile = 0.01,
    clip_negatives = TRUE,
    min_score_spatial = 0.5
) {
  ion_mode  <- match.arg(ion_mode)
  method    <- match.arg(method)
  transform <- match.arg(transform)

  if (is.null(adducts)) adducts <- default_adducts()
  if (!is.list(pkm) || is.null(pkm$mass) || is.null(pkm$intensity)) {
    stop("`pkm` must be a list containing `mass` and `intensity`.", call. = FALSE)
  }

  if (!is.numeric(pkm$mass)) {
    stop("`pkm$mass` must be a numeric vector.", call. = FALSE)
  }

  if (!is.matrix(pkm$intensity)) {
    stop("`pkm$intensity` must be a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(pkm$intensity)) {
    stop("`pkm$intensity` must contain numeric values.", call. = FALSE)
  }

  if (length(pkm$mass) != ncol(pkm$intensity)) {
    stop(
      "Length of `pkm$mass` must match the number of columns in `pkm$intensity`.",
      call. = FALSE
    )
  }

  if (!is.data.frame(adducts)) {
    stop("`adducts` must be a data.frame.", call. = FALSE)
  }

  required_adduct_cols <- c("name", "mode", "mass")
  missing_adduct_cols <- setdiff(required_adduct_cols, names(adducts))

  if (length(missing_adduct_cols) > 0) {
    stop(
      "`adducts` is missing required columns: ",
      paste(missing_adduct_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(adducts$mass)) {
    stop("`adducts$mass` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(tol_ppm) || length(tol_ppm) != 1L ||
      !is.finite(tol_ppm) || tol_ppm <= 0) {
    stop("`tol_ppm` must be a positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(neutral_tol_ppm) || length(neutral_tol_ppm) != 1L ||
      !is.finite(neutral_tol_ppm) || neutral_tol_ppm <= 0) {
    stop("`neutral_tol_ppm` must be a positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(min_quantile) || length(min_quantile) != 1L ||
      !is.finite(min_quantile) || min_quantile < 0 || min_quantile >= 1) {
    stop("`min_quantile` must be a numeric value in the range [0, 1).", call. = FALSE)
  }

  if (!is.logical(clip_negatives) || length(clip_negatives) != 1L) {
    stop("`clip_negatives` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(min_score_spatial) || length(min_score_spatial) != 1L ||
      !is.finite(min_score_spatial) ||
      min_score_spatial < 0 || min_score_spatial > 1) {
    stop("`min_score_spatial` must be a numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.null(feature_idx) &&
      !is.numeric(feature_idx) &&
      !is.integer(feature_idx)) {
    stop(
      "`feature_idx` must be an integer or numeric vector of feature indices.",
      call. = FALSE
    )
  }
  # Gaussian score from a non-negative error
  gaussian_score <- function(err, tol) {
    err <- pmax(0, err)
    exp(-(err^2) / (2 * tol^2))
  }

  empty_adduct_edges <- function() {
    data.frame(
      idx_i = integer(),
      mz_i = numeric(),
      adduct_i = character(),
      idx_j = integer(),
      mz_j = numeric(),
      adduct_j = character(),
      delta_theo = numeric(),
      delta_obs = numeric(),
      delta_err_ppm = numeric(),
      neutral_mass_i = numeric(),
      neutral_mass_j = numeric(),
      neutral_mass_mean = numeric(),
      neutral_err_ppm = numeric(),
      score_spatial = numeric(),
      score_mass = numeric(),
      score_adduct = numeric(),
      is_valid_adduct = logical(),
      stringsAsFactors = FALSE
    )
  }

  # Preprocess two intensity vectors
  preprocess_xy <- function(x, y) {
    keep <- is.finite(x) & is.finite(y) & (x != 0) & (y != 0)
    if (sum(keep) < 3L) return(NULL)

    if (min_quantile > 0) {
      qx <- stats::quantile(x[keep], probs = min_quantile, na.rm = TRUE, type = 7)
      qy <- stats::quantile(y[keep], probs = min_quantile, na.rm = TRUE, type = 7)
      keep <- keep & (x > qx) & (y > qy)
    }

    if (sum(keep) < 3L) return(NULL)

    xk <- x[keep]
    yk <- y[keep]

    if (clip_negatives) {
      xk <- pmax(xk, 0)
      yk <- pmax(yk, 0)
    }

    tf <- switch(
      transform,
      "none"  = identity,
      "log1p" = log1p,
      "zscore" = function(v) {
        sdv <- stats::sd(v)
        if (is.na(sdv) || sdv == 0) rep(0, length(v)) else (v - mean(v)) / sdv
      }
    )

    list(x = tf(xk), y = tf(yk))
  }

  # Similarity score in [0,1]
  score_core <- function(x, y) {
    if (length(x) != length(y) || length(x) < 3L) return(NA_real_)
    if (stats::var(x) == 0 || stats::var(y) == 0) return(NA_real_)

    switch(
      method,
      "pearson" = {
        r <- suppressWarnings(stats::cor(x, y, method = "pearson"))
        if (is.na(r)) NA_real_ else max(0, min(1, r * r))
      },
      "cosine" = {
        num <- sum(x * y)
        den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
        if (den == 0) NA_real_ else max(0, min(1, num / den))
      },
      "spearman" = {
        r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
        if (is.na(r) || r <= 0) 0 else max(0, min(1, r * r))
      }
    )
  }

  mz <- as.numeric(pkm$mass)
  Imat <- pkm$intensity

  # Keep original feature indices, but sort internally by m/z
  ord <- order(mz)
  mz_sorted <- mz[ord]
  Imat_sorted <- Imat[, ord, drop = FALSE]
  original_idx <- ord

  if (is.null(feature_idx)) {
    keep_sorted <- seq_along(mz_sorted)
  } else {
    feature_idx <- as.integer(feature_idx)
    feature_idx <- feature_idx[feature_idx >= 1 & feature_idx <= length(mz)]
    keep_sorted <- which(original_idx %in% feature_idx)
  }

  mz_sorted <- mz_sorted[keep_sorted]
  Imat_sorted <- Imat_sorted[, keep_sorted, drop = FALSE]
  original_idx <- original_idx[keep_sorted]

  p <- length(mz_sorted)

  if (p < 2L) {
    return(empty_adduct_edges())
  }

  # Keep only adducts of the selected ionization mode
  adf <- adducts[adducts$mode == ion_mode, , drop = FALSE]
  adf <- adf[order(adf$mass), , drop = FALSE]
  rownames(adf) <- NULL
  if (nrow(adf) < 2L) {
    stop("At least two adducts are required for ion_mode = '", ion_mode, "'.")
  }

  # Build all unordered adduct pairs
  adduct_pairs <- utils::combn(seq_len(nrow(adf)), 2)

  rows <- list()
  rr <- 0L

  # Helper: find candidate peaks around a target m/z
  get_candidates <- function(target, mz_vec) {
    tol_da <- target * tol_ppm * 1e-6
    which(abs(mz_vec - target) <= tol_da)
  }

  for (k in seq_len(ncol(adduct_pairs))) {
    ia <- adduct_pairs[1, k]
    ib <- adduct_pairs[2, k]

    adduct_a <- adf[ia, , drop = FALSE]
    adduct_b <- adf[ib, , drop = FALSE]

    delta_ab <- adduct_b$mass - adduct_a$mass

    # Search in the forward direction only
    for (i in seq_len(p - 1L)) {
      mz_i <- mz_sorted[i]
      target_j <- mz_i + delta_ab

      cand_j <- get_candidates(target_j, mz_sorted)
      cand_j <- cand_j[cand_j > i]
      if (!length(cand_j)) next

      for (j in cand_j) {
        # spatial/intensity similarity
        pp <- preprocess_xy(Imat_sorted[, i], Imat_sorted[, j])
        if (is.null(pp)) {
          score_spatial <- NA_real_
        } else { score_spatial<- score_core(pp$x, pp$y)}

        mz_j <- mz_sorted[j]
        delta_obs <- mz_j - mz_i
        delta_err_ppm <- 1e6 * abs(delta_obs - delta_ab) / pmax(abs(delta_ab), 1e-12)

        # infer neutral masses using the same convention as default_adducts()
        neutral_i <- mz_i - adduct_a$mass
        neutral_j <- mz_j - adduct_b$mass
        neutral_mean <- mean(c(neutral_i, neutral_j))
        neutral_err_ppm <- 1e6 * abs(neutral_i - neutral_j) / pmax(abs(neutral_mean), 1e-12)

        score_mass <- gaussian_score(neutral_err_ppm, neutral_tol_ppm)

        # Mass is already used as a hard gate through tol_ppm.
        # Therefore, the adduct score is defined as the spatial similarity
        # among mass-compatible adduct pairs.
        score_adduct <- score_spatial

        rr <- rr + 1L
        rows[[rr]] <- data.frame(
          idx_i = original_idx[i],
          mz_i = mz_i,
          adduct_i = adduct_a$name,
          idx_j = original_idx[j],
          mz_j = mz_j,
          adduct_j = adduct_b$name,
          delta_theo = delta_ab,
          delta_obs = delta_obs,
          delta_err_ppm = delta_err_ppm,
          neutral_mass_i = neutral_i,
          neutral_mass_j = neutral_j,
          neutral_mass_mean = neutral_mean,
          neutral_err_ppm = neutral_err_ppm,
          score_spatial = score_spatial,
          score_mass = score_mass,
          score_adduct = score_adduct,
          is_valid_adduct = !is.na(score_adduct) & score_adduct >= min_score_spatial,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(rows)) {
    return(empty_adduct_edges())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  # Keep all compatible adduct hypotheses. Ambiguous feature pairs are resolved
  # downstream during adduct-family grouping.
  out
}
#adduct_res_pos <- adduct_candidates(
#  pkm,
#  ion_mode = "pos"
#)
