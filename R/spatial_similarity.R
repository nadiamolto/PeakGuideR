#' Preprocess two intensity vectors for spatial similarity scoring
#'
#' @description
#' Shared internal preprocessing used by [adduct_candidates()] and
#' [build_single_adduct_candidates()] to prepare a pair of feature intensity
#' vectors before computing a spatial/intensity similarity score.
#'
#' @param x,y Numeric intensity vectors (same length).
#' @param min_quantile Numeric in the range 0 to 1. Feature-wise low-intensity
#'   quantile filter applied after retaining co-detected pixels.
#' @param clip_negatives Logical. If `TRUE`, negative intensities are
#'   truncated to zero before transformation.
#' @param transform Intensity transformation: `"none"`, `"log1p"` or
#'   `"zscore"`.
#'
#' @return A list with `x` and `y` (preprocessed vectors), or `NULL` if fewer
#'   than 3 valid pixels remain.
#' @keywords internal
preprocess_xy <- function(x, y, min_quantile = 0.01, clip_negatives = TRUE,
                           transform = c("none", "log1p", "zscore")) {
  transform <- match.arg(transform)

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

#' Spatial/intensity similarity score in the range 0 to 1
#'
#' @description
#' Shared internal similarity metric used by [adduct_candidates()] and
#' [build_single_adduct_candidates()] on preprocessed intensity vectors (see
#' [preprocess_xy()]).
#'
#' @param x,y Numeric vectors, typically the output of [preprocess_xy()].
#' @param method Similarity metric: `"pearson"`, `"cosine"` or `"spearman"`.
#'
#' @return A single numeric value in the range 0 to 1, or `NA_real_`.
#' @keywords internal
score_core <- function(x, y, method = c("pearson", "cosine", "spearman")) {
  method <- match.arg(method)

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
