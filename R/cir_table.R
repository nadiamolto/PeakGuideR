#' CIR theoretical lookup table (100-1200 Da)
#'
#' Precomputed M+1/M0 ratios for monoisotopic masses from reference formulas.
#'
#' @format ## cir_table
#' A data frame with 1,101 rows and 2 columns:
#' \describe{
#'   \item{mz}{monoisotopic m/z (100-1200 Da, 1 Da steps)}
#'   \item{R_theo}{theoretical CIR ratio M/M+1}
#' }
#'
#' @details
#' Generated via GAM (k=20, gamma=1.2) fitted over carbon counts from
#' ~184k ChEBI/HMDB formulas using natural 13C abundance (p=0.0107).
#'
#' Used by `cir_ratio()` and `cir_score()` for C13 isotope validation.
#'
#' @source ChEBI/HMDB reference databases
"cir_table"
