#' Mass defect of a neutral mass
#'
#' @description
#' Computes the mass defect of neutral masses as the decimal part of the mass,
#' `mass - floor(mass)`, a value in `[0, 1)`. The integer part, `floor(mass)`,
#' is the integer mass the neutral mass belongs to (e.g. 130.0 and 130.9 both
#' belong to 130).
#'
#' @param mass Numeric vector of neutral masses.
#'
#' @return Numeric vector of mass defects (same length as `mass`).
#'
#' @keywords internal
mass_defect_value <- function(mass) {
  mass - floor(mass)
}


#' Reference limits of the mass defect for a neutral mass
#'
#' @description
#' Looks up, in the internal mass defect reference table, the limits of the
#' mass defect observed among HMDB and ChEBI compounds of the same integer
#' mass. The table has one row per integer mass; each neutral mass uses the row
#' of `floor(mass)`.
#'
#' Masses outside the range of the table (55-1300 Da), integer masses without
#' reference compounds and missing masses return `NA` limits.
#'
#' @param mass Numeric vector of neutral masses.
#'
#' @return A data.frame with one row per element of `mass` and columns
#'   `lower` and `upper`.
#'
#' @keywords internal
mass_defect_limits <- function(mass) {
  tab <- mass_defect_table
  integer_mass <- floor(mass)
  row <- integer_mass - min(tab$mass) + 1L

  evaluable <- is.finite(mass) & integer_mass >= min(tab$mass) & integer_mass <= max(tab$mass)
  row[!evaluable] <- NA_integer_

  data.frame(
    lower = tab$lower[row],
    upper = tab$upper[row]
  )
}


#' Mass defect flag for neutral masses
#'
#' @description
#' Flags whether the mass defect of a neutral mass falls within the range of
#' mass defects observed among HMDB and ChEBI compounds of the same integer
#' mass. The flag is informative only: it does not change any score.
#'
#' The mass defect of each neutral mass is computed with
#' [mass_defect_value()] and compared with the limits returned by
#' [mass_defect_limits()].
#'
#' @param neutral_mass Numeric vector of neutral masses.
#'
#' @return A logical vector (same length as `neutral_mass`): `TRUE` when the
#'   mass defect is within the limits, `FALSE` when it is outside, and `NA`
#'   when the mass is missing, outside the range of the reference table
#'   (55-1300 Da) or belongs to an integer mass without reference compounds.
#'
#' @keywords internal
mass_defect_ok <- function(neutral_mass) {
  limits <- mass_defect_limits(neutral_mass)
  defect <- mass_defect_value(neutral_mass)

  defect >= limits$lower & defect <= limits$upper
}


#' Append the mass defect flag to a candidate table
#'
#' @description
#' Adds a `mass_defect_ok` column at the end of a table with a
#' `neutral_mass_consensus` column (`neutral_mass_candidates` or
#' `candidate_annotations`), using [mass_defect_ok()]. The remaining columns
#' and rows are left untouched. If the table has no `neutral_mass_consensus`
#' column, it is returned unchanged.
#'
#' @param x A data.frame with a `neutral_mass_consensus` column.
#'
#' @return `x` with the extra column `mass_defect_ok`.
#'
#' @keywords internal
add_mass_defect_flag <- function(x) {
  if (!is.data.frame(x) || !"neutral_mass_consensus" %in% names(x)) {
    return(x)
  }

  x$mass_defect_ok <- mass_defect_ok(x$neutral_mass_consensus)
  x
}
