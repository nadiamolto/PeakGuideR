
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
#'
#' @examples
#' compound_db <- load_compound_mass_database(quiet = TRUE)
#' head(compound_db)
#'
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


#' Resolve compound identity within a compound mass database
#'
#' @description
#' Assigns a deterministic `compound_identity_id` to each row of `compound_db`
#' by grouping rows that represent the same underlying compound, so that
#' downstream matching (for example against `standards_db`) is never performed
#' against duplicate rows for the same compound.
#'
#' Rows are grouped by exact `StdInChIKey` match first. Rows without a usable
#' `StdInChIKey` are then grouped by exact `StdInChI` match. Rows lacking both
#' identifiers each receive their own identity, since there is no reliable
#' identifier to group them on (molecular formula alone is not treated as
#' proof of identity; see [resolve_cross_source_identity()]).
#'
#' This does **not** use the existing `dup` column in `compound_db`, whose
#' derivation criteria are unknown and unverifiable; deduplication is always
#' recomputed deterministically from `StdInChIKey`/`StdInChI`.
#'
#' @param compound_db A compound mass database data.frame. If it contains
#'   `StdInChIKey` and/or `StdInChI` columns, these are used for grouping.
#'
#' @return `compound_db` with an added integer column `compound_identity_id`.
#'
#' @examples
#' compound_db <- load_compound_mass_database(quiet = TRUE)
#' compound_db <- resolve_compound_identity(compound_db)
#' length(unique(compound_db$compound_identity_id))
#'
#' @export
resolve_compound_identity <- function(compound_db) {
  stopifnot(is.data.frame(compound_db))

  n <- nrow(compound_db)

  clean_id <- function(x) {
    if (is.null(x)) return(rep(NA_character_, n))
    x <- as.character(x)
    x[!nzchar(x)] <- NA_character_
    x
  }

  inchikey <- clean_id(compound_db[["StdInChIKey"]])
  inchi <- clean_id(compound_db[["StdInChI"]])

  key <- ifelse(
    !is.na(inchikey),
    paste0("IK:", inchikey),
    ifelse(
      !is.na(inchi),
      paste0("ICHI:", inchi),
      paste0("ROW:", seq_len(n))
    )
  )

  compound_db$compound_identity_id <- match(key, unique(key))
  compound_db
}


#' Deduplicate a compound mass database by resolved compound identity
#'
#' @description
#' Collapses `compound_db` to one representative row per
#' `compound_identity_id` (see [resolve_compound_identity()]), preferring the
#' row with the most complete identifiers when several rows share the same
#' identity (for example the same compound reported under `Source == "HMDB"`
#' and `Source == "NORMAN"`).
#'
#' @param compound_db A compound mass database data.frame. If it does not
#'   already contain a `compound_identity_id` column,
#'   [resolve_compound_identity()] is applied first.
#'
#' @return A deduplicated data.frame with one row per `compound_identity_id`.
#' @keywords internal
deduplicate_compound_identity <- function(compound_db) {
  stopifnot(is.data.frame(compound_db))

  if (!"compound_identity_id" %in% names(compound_db)) {
    compound_db <- resolve_compound_identity(compound_db)
  }

  has_val <- function(x) {
    if (is.null(x)) return(rep(FALSE, nrow(compound_db)))
    x <- as.character(x)
    !is.na(x) & nzchar(x)
  }

  completeness <- has_val(compound_db[["StdInChIKey"]]) +
    has_val(compound_db[["StdInChI"]]) +
    has_val(compound_db[["SMILES"]])

  compound_db$.completeness <- completeness
  compound_db$.row_order <- seq_len(nrow(compound_db))

  compound_db <- compound_db |>
    dplyr::arrange(
      compound_identity_id,
      dplyr::desc(.completeness),
      .row_order
    ) |>
    dplyr::distinct(compound_identity_id, .keep_all = TRUE) |>
    dplyr::arrange(.row_order)

  compound_db[, setdiff(names(compound_db), c(".completeness", ".row_order")), drop = FALSE]
}


#' Inspect identifier coverage of a compound mass database
#'
#' @description
#' Reports, per `Source`, the percentage of rows with a non-empty
#' `StdInChIKey`, `StdInChI` and `SMILES`. This is useful before relying on
#' identifier-coverage assumptions (for example that a given `Source` never
#' has `StdInChIKey`), since coverage may differ between the bundled example
#' database and the full non-commercial database.
#'
#' @param compound_db A compound mass database data.frame.
#'
#' @return A data.frame with one row per `Source`, columns `n` (row count) and
#'   `pct_StdInChIKey`, `pct_StdInChI`, `pct_SMILES` (percentage of non-empty
#'   values, rounded to 1 decimal).
#'
#' @examples
#' compound_db <- load_compound_mass_database(quiet = TRUE)
#' inspect_database_coverage(compound_db)
#'
#' @export
inspect_database_coverage <- function(compound_db) {
  stopifnot(is.data.frame(compound_db))

  id_cols <- intersect(
    c("StdInChIKey", "StdInChI", "SMILES"),
    names(compound_db)
  )

  source_vec <- if ("Source" %in% names(compound_db)) {
    as.character(compound_db$Source)
  } else {
    rep(NA_character_, nrow(compound_db))
  }

  pct_nonempty <- function(x) {
    x <- as.character(x)
    round(100 * mean(!is.na(x) & nzchar(x)), 1)
  }

  df <- compound_db
  df$.source_group <- ifelse(is.na(source_vec), "NA", source_vec)

  out <- df |>
    dplyr::group_by(Source = .source_group) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(
        dplyr::all_of(id_cols),
        pct_nonempty,
        .names = "pct_{.col}"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(Source)

  for (missing_col in setdiff(paste0("pct_", c("StdInChIKey", "StdInChI", "SMILES")), names(out))) {
    out[[missing_col]] <- NA_real_
  }

  out
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
#'
#' @examples
#' standards_db <- load_standards_adduct_library(quiet = TRUE)
#' head(standards_db)
#'
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
