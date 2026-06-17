#' Mass deviation scoring for isotopic candidate pairs
#'
#' @description
#' Adds a continuous mass-accuracy score to isotopic pairs detected by
#' \code{iso_morphology_candidates()}. This does NOT search new candidates:
#' it only scores the already-detected pairs using their mass error
#' (\code{mass_err_da} / \code{mass_err_ppm}).
#'
#' This is useful when the morphology step uses a hard tolerance window
#' (pass/fail), but you still want a smooth penalty to rank candidates by
#' how close they are to the theoretical delta.
#'
#' @details
#' The function expects \code{pairs} to contain \code{mass_err_ppm} (preferred)
#' or \code{mass_err_da}. The score is computed as either:
#' \itemize{
#'   \item \strong{Gaussian (default)}:
#'     \deqn{score = exp(-(err^2)/(2*sigma^2))}
#'     where \code{sigma = tol/2}.
#'   \item \strong{Linear}:
#'     \deqn{score = pmax(0, 1 - |err|/tol)}
#' }
#'
#' Optionally, a combined score can be produced:
#' \deqn{score\_combined = w\_morph * score\_final + (1-w\_morph) * mass\_dev\_score}
#'
#' @param pairs Data.frame output of \code{iso_morphology_candidates()}.
#'   Must contain \code{score_final} and either \code{mass_err_ppm} or \code{mass_err_da}.
#' @param tol_ppm Numeric. PPM tolerance used as the reference window. Default \code{5}.
#'   (Typically the same \code{tol_ppm} used in \code{iso_morphology_candidates}.)
#' @param tol_da Numeric. Da tolerance used if scoring in Da (only used if \code{use="da"}).
#' @param use Character: \code{"ppm"} (default) or \code{"da"} to decide which error column to use.
#' @param kernel Character: \code{"gaussian"} (default) or \code{"linear"}.
#' @param combine Logical. If \code{TRUE}, compute \code{score_combined}. Default \code{TRUE}.
#' @param w_morph Numeric in range 0 to 1. Weight for morphology when combining. Default \code{0.8}.
#'
#' @return The input \code{pairs} with additional columns:
#' \itemize{
#'   \item \code{mass_dev_err}: absolute mass error in chosen units (ppm or Da)
#'   \item \code{mass_dev_score}: continuous score in the range from 0 to 1.
#'   \item \code{score_combined}: optional combined score
#' }
#'
#' @examples
#' \dontrun{
#' pairs <- iso_morphology_candidates(pm, tol_ppm=5)
#' pairs2 <- iso_mass_deviation_score(pairs, tol_ppm=5, kernel="gaussian", w_morph=0.8)
#' head(pairs2)
#' }
#'
#' @export
iso_mass_deviation_score <- function(pairs,
                                     tol_ppm = 5,
                                     tol_da  = NULL,
                                     use     = c("ppm","da"),
                                     kernel  = c("gaussian","linear"),
                                     combine = TRUE,
                                     w_morph = 0.8) {

  use    <- match.arg(use)
  kernel <- match.arg(kernel)

  stopifnot(is.data.frame(pairs))
  stopifnot("score_final" %in% names(pairs))

  if (use == "ppm") {
    if (!"mass_err_ppm" %in% names(pairs)) stop("pairs must contain 'mass_err_ppm' when use='ppm'.")
    tol <- as.numeric(tol_ppm)
    err <- abs(as.numeric(pairs$mass_err_ppm))
  } else {
    if (!"mass_err_da" %in% names(pairs)) stop("pairs must contain 'mass_err_da' when use='da'.")
    if (is.null(tol_da)) stop("Provide tol_da when use='da'.")
    tol <- as.numeric(tol_da)
    err <- abs(as.numeric(pairs$mass_err_da))
  }

  if (!is.finite(tol) || tol <= 0) stop("Tolerance must be a positive number.")

  mass_dev_score <- if (kernel == "gaussian") {
    sigma <- tol / 2
    exp(-(err^2) / (2 * sigma^2))
  } else {
    pmax(0, 1 - (err / tol))
  }

  out <- pairs
  out$mass_dev_err   <- err
  out$mass_dev_score <- mass_dev_score

  if (isTRUE(combine)) {
    w_morph <- max(0, min(1, as.numeric(w_morph)))
    out$score_combined <- w_morph * out$score_final + (1 - w_morph) * out$mass_dev_score
  }

  out
}
