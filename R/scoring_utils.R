#' Gaussian plausibility score from a non-negative error
#'
#' This script is designed to make the scores compatible.
#'
#' @param err Numeric vector of non-negative errors.
#' @param tol Numeric tolerance.controlling the width of the Gaussian penalty.

#'With this formulation, the score is approximately 0.61 when `err = tol`.
#'
#' @return Numeric vector in range 0 to 1.
#' @keywords internal
gaussian_score <- function(err, tol) {
  err <- pmax(0, err)
  exp(-(err^2) / (2 * tol^2))
}

#' Convert a ppm mass error into a 0-1 plausibility score
#'
#' @description
#' Converts a candidate mass error (in ppm) into a `mass_error_score` in the
#' range 0 to 1, using the same Gaussian penalty as [gaussian_score()].
#'
#' @param candidate_ppm_error Numeric vector of non-negative ppm errors.
#' @param ppm_tol Numeric tolerance (ppm) controlling the width of the penalty.
#'
#' @return Numeric vector in range 0 to 1.
#' @keywords internal
ppm_error_to_score <- function(candidate_ppm_error, ppm_tol) {
  gaussian_score(candidate_ppm_error, ppm_tol)
}
