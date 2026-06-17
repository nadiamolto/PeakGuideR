
#' Load the PeakGuideR compound mass database
#'
#' @description
#' Loads the example compound mass database included with PeakGuideR.
#' The full non-commercial database is distributed separately and can be
#' supplied manually to the workflow through the `compound_db` argument.
#'
#' @param quiet Logical. If `FALSE`, prints a message.
#'
#' @return A data.frame.
#' @export
load_compound_mass_database <- function(quiet = FALSE) {
  path <- system.file(
    "extdata",
    "compound_mass_database_example.rds",
    package = "PeakGuideR"
  )

  if (!nzchar(path)) {
    stop(
      "Could not find 'compound_mass_database_example.rds' in package extdata.",
      call. = FALSE
    )
  }

  if (!isTRUE(quiet)) {
    message(
      "Loading the example compound mass database included with PeakGuideR."
    )
  }

  readRDS(path)
}


#' Load the PeakGuideR standard adduct library
#'
#' @description
#' Loads the example standard adduct library included with PeakGuideR.
#' The full non-commercial standard adduct library is distributed separately
#' and can be supplied manually to the workflow through the `standards_db`
#' argument.
#'
#' @param quiet Logical. If `FALSE`, prints a message.
#'
#' @return A data.frame.
#' @export
load_standards_adduct_library <- function(quiet = FALSE) {
  path <- system.file(
    "extdata",
    "standards_adduct_library_example.rds",
    package = "PeakGuideR"
  )

  if (!nzchar(path)) {
    stop(
      "Could not find 'standards_adduct_library_example.rds' in package extdata.",
      call. = FALSE
    )
  }

  if (!isTRUE(quiet)) {
    message(
      "Loading the example standard adduct library included with PeakGuideR."
    )
  }

  readRDS(path)
}
