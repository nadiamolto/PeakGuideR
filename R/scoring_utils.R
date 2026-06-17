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
