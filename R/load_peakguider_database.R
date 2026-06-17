#' Load external PeakGuideR databases
#'
#' @description
#' Loads full PeakGuideR annotation databases downloaded separately, for example
#' from Zenodo.
#'
#' @param compound_db_path Path to `compound_mass_database_noncommercial.rds`.
#' @param standards_db_path Optional path to
#'   `standards_adduct_library_noncommercial.rds`.
#'
#' @return A list with `compound_db` and `standards_db`.
#' @export
load_peakguider_databases <- function(
    compound_db_path,
    standards_db_path = NULL
) {
  if (missing(compound_db_path) || !nzchar(compound_db_path)) {
    stop("Please provide `compound_db_path`.", call. = FALSE)
  }

  if (!file.exists(compound_db_path)) {
    stop("Could not find compound database file: ", compound_db_path, call. = FALSE)
  }

  compound_db <- readRDS(compound_db_path)

  if (!is.data.frame(compound_db)) {
    stop("The compound database must be a data.frame.", call. = FALSE)
  }

  standards_db <- NULL

  if (!is.null(standards_db_path)) {
    if (!file.exists(standards_db_path)) {
      stop("Could not find standards database file: ", standards_db_path, call. = FALSE)
    }

    standards_db <- readRDS(standards_db_path)

    if (!is.data.frame(standards_db)) {
      stop("The standards database must be a data.frame.", call. = FALSE)
    }
  }

  list(
    compound_db = compound_db,
    standards_db = standards_db
  )
}
