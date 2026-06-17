#' Theoretical CIR ratio
#'
#' @param mz monoisotopic masses (numeric vector).
#' @return CIR theoretical ratios (same length).
#' @export
cir_ratio <- function(mz) {
  stats::approx(cir_table$mz, cir_table$R_theo, xout = mz, rule = 2)$y
}
