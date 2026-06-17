#' EIPS theoretical lookup table
#'
#' Precomputed expected isotope presence scores (EIPS) for different isotope
#' types and atom counts.
#'
#' @format A data.frame with columns used to store theoretical or empirical
#' isotope presence score expectations.
#'
#' @source Internal PeakGuideR reference table.
#'
#' @keywords datasets
#' @name eips_table
#' @docType data
#' @usage data(eips_table)
NULL


#' EIPS elemental isotope-pattern reference table
#'
#' Reference table used by `eips_score()` to evaluate elemental isotope-pattern
#' support. The EIPS scoring step can use isotope-pattern evidence from several
#' elements, including Br, Cl, N, O and S.
#'
#' @format A data.frame with precomputed elemental isotope-pattern reference
#' values used by PeakGuideR.
#'
#' @source Internal PeakGuideR reference table.
#'
#' @keywords datasets
#' @name eips_n_table
#' @docType data
#' @usage data(eips_n_table)
NULL
