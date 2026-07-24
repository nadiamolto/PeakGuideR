#' Cluster neutral masses by ppm tolerance
#'
#' @description
#' Shared internal sequential clustering used to assign `neutral_mass_id`
#' groups from a numeric vector of inferred neutral masses. Masses are
#' processed in ascending order; a new cluster starts whenever the ppm
#' distance to the current cluster's anchor mass (its first, lowest member)
#' exceeds `neutral_cluster_ppm`.
#'
#' Used by both [build_neutral_mass_candidates()] and
#' [build_single_adduct_candidates()] so that, given the same adduct-family
#' input and `neutral_cluster_ppm`, they assign identical `neutral_mass_id`
#' values to family-derived neutral masses.
#'
#' @param masses Numeric vector of neutral masses.
#' @param neutral_cluster_ppm PPM tolerance.
#'
#' @return Integer vector of cluster ids, same length and order as `masses`.
#' @keywords internal
assign_neutral_mass_clusters <- function(masses, neutral_cluster_ppm) {
  n <- length(masses)
  cluster_id <- integer(n)
  if (n == 0L) return(cluster_id)

  ord <- order(masses)
  sorted <- masses[ord]

  cluster_id_sorted <- integer(n)
  current_id <- 1L
  cluster_id_sorted[1] <- current_id
  current_mass <- sorted[1]

  if (n > 1L) {
    for (i in 2:n) {
      ppm_diff <- 1e6 * abs(sorted[i] - current_mass) /
        pmax(abs(current_mass), 1e-12)

      if (is.finite(ppm_diff) && ppm_diff <= neutral_cluster_ppm) {
        cluster_id_sorted[i] <- current_id
      } else {
        current_id <- current_id + 1L
        cluster_id_sorted[i] <- current_id
        current_mass <- sorted[i]
      }
    }
  }

  cluster_id[ord] <- cluster_id_sorted
  cluster_id
}


#' Union detected adducts per neutral mass id
#'
#' @description
#' Combines one or more data.frames (each with at least `neutral_mass_id` and
#' `adduct` columns) into a single lookup of which adducts were detected
#' computationally for each `neutral_mass_id`. Used to build the
#' `detected_adducts_set` consumed by [compute_standard_adduct_recovery()].
#'
#' @param ... Data.frames with `neutral_mass_id` and `adduct` columns, or
#'   `NULL` (ignored).
#'
#' @return A named list keyed by `as.character(neutral_mass_id)`, each element
#'   a character vector of unique detected adducts.
#' @keywords internal
compute_detected_adducts_map <- function(...) {
  dfs <- list(...)
  dfs <- Filter(function(d) is.data.frame(d) && nrow(d) > 0, dfs)

  if (!length(dfs)) return(list())

  dfs <- lapply(dfs, function(d) {
    data.frame(
      neutral_mass_id = as.character(d[["neutral_mass_id"]]),
      adduct = as.character(d[["adduct"]]),
      stringsAsFactors = FALSE
    )
  })

  combined <- do.call(rbind, dfs)
  combined <- combined[!is.na(combined$neutral_mass_id) & !is.na(combined$adduct), , drop = FALSE]

  if (!nrow(combined)) return(list())

  split(combined$adduct, combined$neutral_mass_id)
}


#' Match inferred neutral masses against the standard-adduct library by mass
#'
#' @description
#' Retrieves standard-adduct library compounds within `ppm_tol` of each
#' inferred neutral mass, symmetrically to how `compound_db` is matched by
#' mass elsewhere in PeakGuideR.
#'
#' Unlike the previous implementation, this function never joins on
#' identifiers hanging off a `compound_db` candidate row: `standards_db` is
#' treated as an independent evidence source that is retrieved purely by
#' neutral mass, ion mode and (optionally) matrix. Cross-source identity is
#' resolved afterwards by [resolve_cross_source_identity()].
#'
#' Because `standards_db` typically has several adduct rows per compound
#' (all sharing the same `NEUTRAL_MONOISOTOPIC_MASS`), matches are collapsed
#' to one row per standard compound and per `neutral_mass_id`, with
#' `standard_db_adducts` listing every annotated adduct for that compound
#' (used downstream to compute `standard_adduct_recovery_score`). When two or
#' more standard-library compounds (isomers) fall within tolerance of the
#' same neutral mass, they are kept as separate rows rather than merged.
#'
#' @param neutral_masses A data.frame with at least `neutral_mass_id` and
#'   `neutral_mass_consensus`.
#' @param standards_db Standard-adduct library data.frame. Must contain
#'   `NEUTRAL_MONOISOTOPIC_MASS` and `POLARITY`.
#' @param ion_mode `"pos"` or `"neg"`.
#' @param matrix Matrix name. `standards_db` is matrix-specific (the bundled
#'   library only covers HCCA/DEA), so a match against it is only meaningful
#'   when the caller's matrix is known. When supplied and `standards_db`
#'   contains a `matrix` column, only rows whose `matrix` value contains
#'   `matrix` (case-insensitive) are kept. Use `NULL` to skip standard-adduct
#'   library matching entirely (no rows are returned and `standards_db` is
#'   not consulted), rather than comparing against every matrix indiscriminately.
#' @param ppm_tol PPM tolerance for neutral-mass matching.
#' @param quiet Logical. If `FALSE` (the default), prints a message when
#'   `matrix = NULL` explaining that standard-adduct library matching was
#'   skipped.
#'
#' @return A data.frame with one row per `neutral_mass_id` and standard
#'   compound candidate, with columns `neutral_mass_id`, `candidate_source`
#'   (`"standards_db"`), `candidate_name`, `candidate_formula`,
#'   `candidate_neutral_mass`, `candidate_ppm_error`, `candidate_db_id`,
#'   `candidate_inchikey`, `candidate_smiles`, `candidate_hmdb`,
#'   `candidate_chebi` and `standard_db_adducts` (semicolon-separated adducts
#'   annotated for that compound in `standards_db`, after ion-mode/matrix
#'   filtering).
#'
#' @examples
#' standards_db <- load_standards_adduct_library(quiet = TRUE)
#'
#' neutral_masses <- data.frame(
#'   neutral_mass_id = 1L,
#'   neutral_mass_consensus = standards_db$NEUTRAL_MONOISOTOPIC_MASS[1]
#' )
#'
#' match_standards_by_mass(
#'   neutral_masses,
#'   standards_db,
#'   ion_mode = "neg",
#'   matrix = "HCCA"
#' )
#'
#' @export
match_standards_by_mass <- function(
    neutral_masses,
    standards_db,
    ion_mode = c("pos", "neg"),
    matrix = NULL,
    ppm_tol = 5,
    quiet = FALSE
) {
  ion_mode <- match.arg(ion_mode)

  stopifnot(is.data.frame(neutral_masses))
  required_neutral_cols <- c("neutral_mass_id", "neutral_mass_consensus")
  if (!all(required_neutral_cols %in% names(neutral_masses))) {
    stop(
      "`neutral_masses` must contain `neutral_mass_id` and `neutral_mass_consensus`.",
      call. = FALSE
    )
  }

  empty_out <- function() {
    data.frame(
      neutral_mass_id = integer(),
      candidate_source = character(),
      candidate_name = character(),
      candidate_formula = character(),
      candidate_neutral_mass = numeric(),
      candidate_ppm_error = numeric(),
      candidate_db_id = character(),
      candidate_inchikey = character(),
      candidate_smiles = character(),
      candidate_hmdb = character(),
      candidate_chebi = character(),
      standard_db_adducts = character(),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(matrix)) {
    if (!isTRUE(quiet)) {
      message("No matrix specified: skipping standard-adduct library matching.")
    }
    return(empty_out())
  }

  stopifnot(is.data.frame(standards_db))
  required_std_cols <- c("NEUTRAL_MONOISOTOPIC_MASS", "POLARITY")
  missing_std_cols <- setdiff(required_std_cols, names(standards_db))
  if (length(missing_std_cols) > 0) {
    stop(
      "`standards_db` is missing required columns: ",
      paste(missing_std_cols, collapse = ", "),
      call. = FALSE
    )
  }

  std <- standards_db |>
    dplyr::filter(POLARITY == ion_mode)

  if ("matrix" %in% names(std)) {
    std <- std[grepl(matrix, std[["matrix"]], ignore.case = TRUE), , drop = FALSE]
  }

  if (nrow(std) == 0) {
    return(empty_out())
  }

  optional_cols <- c(
    "COMPOUND_ID", "Master_List_NAME", "HMDB_clean", "ChEBI",
    "InCHIKey", "SMILES", "MOLECULAR_FORMULA", "adduct"
  )
  for (col in optional_cols) {
    if (!col %in% names(std)) std[[col]] <- NA_character_
  }

  std2 <- std |>
    dplyr::mutate(
      .std_compound_id = as.character(COMPOUND_ID),
      .std_name = as.character(Master_List_NAME),
      .std_hmdb = as.character(HMDB_clean),
      .std_chebi = as.character(ChEBI),
      .std_inchikey = as.character(InCHIKey),
      .std_smiles = as.character(SMILES),
      .std_formula = as.character(MOLECULAR_FORMULA),
      .std_adduct = as.character(adduct),
      candidate_neutral_mass = as.numeric(NEUTRAL_MONOISOTOPIC_MASS)
    ) |>
    dplyr::filter(is.finite(candidate_neutral_mass))

  if (nrow(std2) == 0) {
    return(empty_out())
  }

  neutral_for_matching <- neutral_masses |>
    dplyr::select(neutral_mass_id, neutral_mass_consensus) |>
    dplyr::distinct()

  matched <- tidyr::crossing(
    neutral_for_matching,
    std2 |>
      dplyr::select(
        .std_compound_id, .std_name, .std_hmdb, .std_chebi,
        .std_inchikey, .std_smiles, .std_formula, .std_adduct,
        candidate_neutral_mass
      )
  ) |>
    dplyr::mutate(
      candidate_ppm_error =
        1e6 * abs(candidate_neutral_mass - neutral_mass_consensus) /
        pmax(abs(neutral_mass_consensus), 1e-12)
    ) |>
    dplyr::filter(
      is.finite(candidate_ppm_error),
      candidate_ppm_error <= ppm_tol
    )

  if (nrow(matched) == 0) {
    return(empty_out())
  }

  std_key <- dplyr::coalesce(matched$.std_compound_id, matched$.std_name)

  matched |>
    dplyr::mutate(.std_key = std_key) |>
    dplyr::group_by(
      neutral_mass_id, .std_key,
      .std_name, .std_formula, .std_inchikey, .std_smiles,
      .std_hmdb, .std_chebi, candidate_neutral_mass
    ) |>
    dplyr::summarise(
      candidate_ppm_error = min(candidate_ppm_error, na.rm = TRUE),
      standard_db_adducts = paste(sort(unique(.std_adduct)), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      neutral_mass_id,
      candidate_source = "standards_db",
      candidate_name = .std_name,
      candidate_formula = .std_formula,
      candidate_neutral_mass,
      candidate_ppm_error,
      candidate_db_id = .std_key,
      candidate_inchikey = .std_inchikey,
      candidate_smiles = .std_smiles,
      candidate_hmdb = .std_hmdb,
      candidate_chebi = .std_chebi,
      standard_db_adducts
    ) |>
    dplyr::arrange(neutral_mass_id, candidate_ppm_error)
}


#' Match inferred neutral masses against a compound mass database
#'
#' @description
#' Retrieves `compound_db` candidates within `ppm_tol` of each inferred
#' neutral mass. `compound_db` is deduplicated by resolved compound identity
#' (see [resolve_compound_identity()] and [deduplicate_compound_identity()])
#' before matching, so that compounds reported under more than one `Source`
#' (for example the same compound in both `HMDB` and `NORMAN`) are matched
#' only once per neutral mass.
#'
#' @param neutral_masses A data.frame with at least `neutral_mass_id` and
#'   `neutral_mass_consensus`.
#' @param compound_db Compound mass database. Must contain a numeric
#'   `MonoisotopicMass` column.
#' @param ppm_tol PPM tolerance for neutral-mass matching.
#' @param top_n Maximum number of candidates kept per `neutral_mass_id`,
#'   ranked by ascending ppm error. Use `NULL` to retain all candidates
#'   within `ppm_tol`.
#'
#' @return A data.frame with one row per `neutral_mass_id` and compound
#'   candidate: `neutral_mass_id`, `candidate_source`, `candidate_db_id`,
#'   `candidate_name`, `candidate_formula`, `candidate_neutral_mass`,
#'   `candidate_inchi`, `candidate_inchikey`, `candidate_smiles`,
#'   `candidate_kegg`, `candidate_mass`, `candidate_ppm_error`,
#'   `compound_identity_id`.
#' @keywords internal
match_compound_db_by_mass <- function(
    neutral_masses,
    compound_db,
    ppm_tol = 5,
    top_n = 10L
) {
  stopifnot(is.data.frame(neutral_masses))
  stopifnot(is.data.frame(compound_db))

  if (!"MonoisotopicMass" %in% names(compound_db)) {
    stop("`compound_db` must contain a `MonoisotopicMass` column.", call. = FALSE)
  }

  optional_cols <- c(
    "Source", "DB_ID", "Name", "MolecularFormula",
    "StdInChI", "StdInChIKey", "SMILES", "Kegg"
  )
  for (col in optional_cols) {
    if (!col %in% names(compound_db)) compound_db[[col]] <- NA_character_
  }

  compound_db <- resolve_compound_identity(compound_db)
  compound_db <- deduplicate_compound_identity(compound_db)

  db2 <- compound_db |>
    dplyr::select(
      dplyr::all_of(c(optional_cols, "MonoisotopicMass", "compound_identity_id"))
    ) |>
    dplyr::mutate(candidate_mass = as.numeric(MonoisotopicMass)) |>
    dplyr::filter(is.finite(candidate_mass))

  neutral_for_matching <- neutral_masses |>
    dplyr::select(neutral_mass_id, neutral_mass_consensus) |>
    dplyr::distinct()

  # Sorted-mass interval search (O(n_neutral_masses * log(n_compounds)))
  # instead of a full cross join, since build_single_adduct_candidates() can
  # generate many distinct neutral masses to match against a large
  # compound_db.
  db2 <- db2[order(db2$candidate_mass), , drop = FALSE]
  database_masses <- db2$candidate_mass

  matched_rows <- vector("list", nrow(neutral_for_matching))

  for (i in seq_len(nrow(neutral_for_matching))) {
    nm_id <- neutral_for_matching$neutral_mass_id[[i]]
    neutral_mass <- neutral_for_matching$neutral_mass_consensus[[i]]

    if (!is.finite(neutral_mass) || !length(database_masses)) next

    absolute_tolerance <- abs(neutral_mass) * ppm_tol * 1e-6
    lower_mass <- neutral_mass - absolute_tolerance
    upper_mass <- neutral_mass + absolute_tolerance

    first_index <- findInterval(lower_mass, database_masses) + 1L
    last_index <- findInterval(upper_mass, database_masses)

    if (first_index > last_index ||
        first_index > length(database_masses) || last_index < 1L) {
      next
    }

    first_index <- max(first_index, 1L)
    last_index <- min(last_index, length(database_masses))

    subset_db <- db2[seq.int(first_index, last_index), , drop = FALSE]

    candidate_ppm_error <- 1e6 * abs(subset_db$candidate_mass - neutral_mass) /
      pmax(abs(neutral_mass), 1e-12)

    keep <- is.finite(candidate_ppm_error) & candidate_ppm_error <= ppm_tol
    if (!any(keep)) next

    subset_db <- subset_db[keep, , drop = FALSE]
    candidate_ppm_error <- candidate_ppm_error[keep]

    matched_rows[[i]] <- dplyr::bind_cols(
      neutral_mass_id = nm_id,
      neutral_mass_consensus = neutral_mass,
      subset_db,
      candidate_ppm_error = candidate_ppm_error
    )
  }

  has_matches <- any(!vapply(matched_rows, is.null, logical(1)))

  out <- if (has_matches) {
    dplyr::bind_rows(matched_rows) |>
      dplyr::arrange(neutral_mass_id, candidate_ppm_error)
  } else {
    dplyr::bind_cols(
      neutral_mass_id = integer(),
      neutral_mass_consensus = numeric(),
      db2[0, , drop = FALSE],
      candidate_ppm_error = numeric()
    )
  }

  if (!is.null(top_n)) {
    out <- out |>
      dplyr::group_by(neutral_mass_id) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  out |>
    dplyr::rename(
      candidate_source = Source,
      candidate_db_id = DB_ID,
      candidate_name = Name,
      candidate_formula = MolecularFormula,
      candidate_neutral_mass = MonoisotopicMass,
      candidate_inchi = StdInChI,
      candidate_inchikey = StdInChIKey,
      candidate_smiles = SMILES,
      candidate_kegg = Kegg
    )
}


#' Resolve compound identity across a broad compound database and the
#' standard-adduct library
#'
#' @description
#' Within each `neutral_mass_id`, links `compound_db`-derived candidates
#' (`compound_candidates`) with `standards_db`-derived candidates
#' (`standard_candidates`, typically from [match_standards_by_mass()]) using
#' a hierarchy of identifiers, stopping at the first level that produces a
#' match: `InChIKey` -> `InChI` -> `SMILES` -> shared database ID (`HMDB` or
#' `ChEBI`, only when the `compound_db` row's `Source` matches that identifier
#' type). Matching is 1:1: a given candidate is linked to at most one
#' counterpart, using the strongest available identifier level.
#'
#' A match on molecular formula alone (with or without an exact name match)
#' is **never** treated as identity. Such rows are kept as separate
#' `"broad_db_only"` / `"standards_only"` rows and cross-referenced through
#' `possible_link_formula_only` for manual review.
#'
#' `identity_match_type` is bookkeeping only: it records which identifier
#' level produced a `"both"` fusion, but must not be used as a numeric
#' evidence score.
#'
#' @param compound_candidates A data.frame of `compound_db`-derived
#'   candidates with at least `neutral_mass_id`, `candidate_source`,
#'   `candidate_db_id`, `candidate_name`, `candidate_formula`,
#'   `candidate_neutral_mass`, `candidate_ppm_error`, `candidate_inchi`,
#'   `candidate_inchikey`, `candidate_smiles`.
#' @param standard_candidates A data.frame of `standards_db`-derived
#'   candidates, typically the output of [match_standards_by_mass()], or
#'   `NULL`/empty if standard-adduct support is not being used.
#'
#' @return A data.frame with one row per resolved candidate identity and
#'   columns: `neutral_mass_id`, `source` (`"both"`, `"broad_db_only"` or
#'   `"standards_only"`), `identity_match_type`, `possible_link_formula_only`,
#'   `candidate_ppm_error`, `candidate_neutral_mass`, `broad_db_name`,
#'   `broad_db_formula`, `broad_db_inchikey`, `broad_db_inchi`,
#'   `broad_db_smiles`, `broad_db_id`, `broad_db_source`,
#'   `standard_db_name`, `standard_db_formula`, `standard_db_inchikey`,
#'   `standard_db_smiles`, `standard_db_compound_id`, `standard_db_adducts`.
#'
#' @examples
#' # "Citric acid" and "Isocitric acid" share the formula C6H8O7, but only
#' # "Citric acid" shares a ChEBI id with the standard "CITRATE": formula
#' # alone must not fuse the two compounds.
#' compound_candidates <- data.frame(
#'   neutral_mass_id = c(1L, 1L),
#'   candidate_source = c("CHEBI", "CHEBI"),
#'   candidate_db_id = c("30769", "30887"),
#'   candidate_name = c("Citric acid", "Isocitric acid"),
#'   candidate_formula = c("C6H8O7", "C6H8O7"),
#'   candidate_neutral_mass = c(192.027, 192.027),
#'   candidate_ppm_error = c(0, 0),
#'   candidate_inchi = NA_character_,
#'   candidate_inchikey = NA_character_,
#'   candidate_smiles = NA_character_,
#'   stringsAsFactors = FALSE
#' )
#'
#' standard_candidates <- data.frame(
#'   neutral_mass_id = 1L,
#'   candidate_source = "standards_db",
#'   candidate_name = "CITRATE",
#'   candidate_formula = "C6H8O7",
#'   candidate_neutral_mass = 192.027,
#'   candidate_ppm_error = 0,
#'   candidate_db_id = "P1A5",
#'   candidate_inchikey = "KRKNYBCHXYNGOX-UHFFFAOYSA-N",
#'   candidate_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O",
#'   candidate_hmdb = "HMDB0000094",
#'   candidate_chebi = "30769",
#'   standard_db_adducts = "[M-H]-",
#'   stringsAsFactors = FALSE
#' )
#'
#' resolve_cross_source_identity(compound_candidates, standard_candidates)
#'
#' @export
resolve_cross_source_identity <- function(
    compound_candidates,
    standard_candidates = NULL
) {
  stopifnot(is.data.frame(compound_candidates))

  norm_chr <- function(x, n = length(x)) {
    if (is.null(x)) return(rep(NA_character_, n))
    x <- as.character(x)
    x[!nzchar(x)] <- NA_character_
    x
  }

  n_cc <- nrow(compound_candidates)

  cc <- compound_candidates
  cc$.broad_id <- seq_len(n_cc)
  cc$.broad_inchikey <- norm_chr(cc[["candidate_inchikey"]], n_cc)
  cc$.broad_inchi <- norm_chr(cc[["candidate_inchi"]], n_cc)
  cc$.broad_smiles <- norm_chr(cc[["candidate_smiles"]], n_cc)
  cc$.broad_formula <- norm_chr(cc[["candidate_formula"]], n_cc)
  cc$.broad_shared_id <- dplyr::case_when(
    toupper(cc[["candidate_source"]]) == "HMDB" ~ norm_chr(cc[["candidate_db_id"]], n_cc),
    toupper(cc[["candidate_source"]]) == "CHEBI" ~ norm_chr(cc[["candidate_db_id"]], n_cc),
    TRUE ~ NA_character_
  )
  cc$.broad_shared_id_type <- dplyr::case_when(
    toupper(cc[["candidate_source"]]) == "HMDB" ~ "hmdb",
    toupper(cc[["candidate_source"]]) == "CHEBI" ~ "chebi",
    TRUE ~ NA_character_
  )

  has_std <- !is.null(standard_candidates) &&
    is.data.frame(standard_candidates) &&
    nrow(standard_candidates) > 0

  broad_only_rows <- function(rows) {
    rows |>
      dplyr::transmute(
        neutral_mass_id,
        source = "broad_db_only",
        identity_match_type = NA_character_,
        possible_link_formula_only = NA_character_,
        candidate_ppm_error,
        candidate_neutral_mass,
        broad_db_name = candidate_name,
        broad_db_formula = candidate_formula,
        broad_db_inchikey = candidate_inchikey,
        broad_db_inchi = candidate_inchi,
        broad_db_smiles = candidate_smiles,
        broad_db_id = candidate_db_id,
        broad_db_source = candidate_source,
        standard_db_name = NA_character_,
        standard_db_formula = NA_character_,
        standard_db_inchikey = NA_character_,
        standard_db_smiles = NA_character_,
        standard_db_compound_id = NA_character_,
        standard_db_adducts = NA_character_
      )
  }

  if (!has_std) {
    return(broad_only_rows(cc))
  }

  n_sc <- nrow(standard_candidates)
  sc <- standard_candidates
  sc$.std_id <- seq_len(n_sc)
  sc$.std_inchikey <- norm_chr(sc[["candidate_inchikey"]], n_sc)
  sc$.std_smiles <- norm_chr(sc[["candidate_smiles"]], n_sc)
  sc$.std_formula <- norm_chr(sc[["candidate_formula"]], n_sc)
  sc$.std_hmdb <- norm_chr(sc[["candidate_hmdb"]], n_sc)
  sc$.std_chebi <- norm_chr(sc[["candidate_chebi"]], n_sc)
  # standards_db carries no InChI column, so the "inchi" hierarchy level
  # never has a counterpart key on this side.
  sc$.std_inchi <- rep(NA_character_, n_sc)

  greedy_pairs <- function(cc_pool, sc_pool, cc_key_col, sc_key_col) {
    cc_key <- cc_pool[[cc_key_col]]
    sc_key <- sc_pool[[sc_key_col]]

    cc_df <- data.frame(
      broad_row = cc_pool$.broad_id,
      neutral_mass_id = cc_pool$neutral_mass_id,
      key = cc_key,
      stringsAsFactors = FALSE
    )
    sc_df <- data.frame(
      std_row = sc_pool$.std_id,
      neutral_mass_id = sc_pool$neutral_mass_id,
      key = sc_key,
      stringsAsFactors = FALSE
    )

    cc_df <- cc_df[!is.na(cc_df$key), , drop = FALSE]
    sc_df <- sc_df[!is.na(sc_df$key), , drop = FALSE]

    if (!nrow(cc_df) || !nrow(sc_df)) {
      return(data.frame(broad_row = integer(), std_row = integer()))
    }

    cand <- merge(cc_df, sc_df, by = c("neutral_mass_id", "key"))
    if (!nrow(cand)) {
      return(data.frame(broad_row = integer(), std_row = integer()))
    }

    cand <- cand[order(cand$broad_row, cand$std_row), , drop = FALSE]

    used_broad <- integer(0)
    used_std <- integer(0)
    keep <- logical(nrow(cand))

    for (i in seq_len(nrow(cand))) {
      b <- cand$broad_row[i]
      s <- cand$std_row[i]
      if (!(b %in% used_broad) && !(s %in% used_std)) {
        keep[i] <- TRUE
        used_broad <- c(used_broad, b)
        used_std <- c(used_std, s)
      }
    }

    cand[keep, c("broad_row", "std_row"), drop = FALSE]
  }

  remaining_broad <- cc$.broad_id
  remaining_std <- sc$.std_id
  matched_pairs <- list()

  levels <- list(
    list(name = "inchikey", cc_col = ".broad_inchikey", sc_col = ".std_inchikey"),
    list(name = "inchi", cc_col = ".broad_inchi", sc_col = ".std_inchi"),
    list(name = "smiles", cc_col = ".broad_smiles", sc_col = ".std_smiles"),
    list(name = "shared_id", cc_col = ".broad_shared_id", sc_col = ".std_hmdb", type = "hmdb"),
    list(name = "shared_id", cc_col = ".broad_shared_id", sc_col = ".std_chebi", type = "chebi")
  )

  for (lvl in levels) {
    cc_pool <- cc[cc$.broad_id %in% remaining_broad, , drop = FALSE]
    sc_pool <- sc[sc$.std_id %in% remaining_std, , drop = FALSE]

    if (!is.null(lvl$type)) {
      cc_pool <- cc_pool[
        is.na(cc_pool$.broad_shared_id_type) | cc_pool$.broad_shared_id_type == lvl$type,
        , drop = FALSE
      ]
    }

    if (!nrow(cc_pool) || !nrow(sc_pool)) next

    pairs <- greedy_pairs(cc_pool, sc_pool, lvl$cc_col, lvl$sc_col)
    if (!nrow(pairs)) next

    pairs$identity_match_type <- lvl$name
    matched_pairs[[length(matched_pairs) + 1L]] <- pairs

    remaining_broad <- setdiff(remaining_broad, pairs$broad_row)
    remaining_std <- setdiff(remaining_std, pairs$std_row)
  }

  both_rows <- NULL
  if (length(matched_pairs)) {
    all_pairs <- do.call(rbind, matched_pairs)

    cc_matched <- cc[match(all_pairs$broad_row, cc$.broad_id), , drop = FALSE]
    sc_matched <- sc[match(all_pairs$std_row, sc$.std_id), , drop = FALSE]

    both_rows <- data.frame(
      neutral_mass_id = cc_matched$neutral_mass_id,
      source = "both",
      identity_match_type = all_pairs$identity_match_type,
      possible_link_formula_only = NA_character_,
      candidate_ppm_error = cc_matched$candidate_ppm_error,
      candidate_neutral_mass = cc_matched$candidate_neutral_mass,
      broad_db_name = cc_matched$candidate_name,
      broad_db_formula = cc_matched$candidate_formula,
      broad_db_inchikey = cc_matched$candidate_inchikey,
      broad_db_inchi = cc_matched$candidate_inchi,
      broad_db_smiles = cc_matched$candidate_smiles,
      broad_db_id = cc_matched$candidate_db_id,
      broad_db_source = cc_matched$candidate_source,
      standard_db_name = sc_matched$candidate_name,
      standard_db_formula = sc_matched$candidate_formula,
      standard_db_inchikey = sc_matched$candidate_inchikey,
      standard_db_smiles = sc_matched$candidate_smiles,
      standard_db_compound_id = sc_matched$candidate_db_id,
      standard_db_adducts = sc_matched$standard_db_adducts,
      stringsAsFactors = FALSE
    )
  }

  cc_left <- cc[cc$.broad_id %in% remaining_broad, , drop = FALSE]
  sc_left <- sc[sc$.std_id %in% remaining_std, , drop = FALSE]

  # Cross-reference unmatched rows that only share a molecular formula, for
  # manual review. This never fuses identity.
  formula_link_broad <- rep(NA_character_, nrow(cc_left))
  formula_link_std <- rep(NA_character_, nrow(sc_left))

  if (nrow(cc_left) && nrow(sc_left)) {
    broad_lbl <- paste0(
      dplyr::coalesce(cc_left$candidate_source, "broad_db"), ":",
      dplyr::coalesce(as.character(cc_left$candidate_db_id), cc_left$candidate_name)
    )
    std_lbl <- paste0(
      "standard:",
      dplyr::coalesce(as.character(sc_left$candidate_db_id), sc_left$candidate_name)
    )

    cc_f <- data.frame(
      broad_row = cc_left$.broad_id,
      neutral_mass_id = cc_left$neutral_mass_id,
      formula = cc_left$.broad_formula,
      label = broad_lbl,
      stringsAsFactors = FALSE
    )
    sc_f <- data.frame(
      std_row = sc_left$.std_id,
      neutral_mass_id = sc_left$neutral_mass_id,
      formula = sc_left$.std_formula,
      label = std_lbl,
      stringsAsFactors = FALSE
    )
    cc_f <- cc_f[!is.na(cc_f$formula), , drop = FALSE]
    sc_f <- sc_f[!is.na(sc_f$formula), , drop = FALSE]

    if (nrow(cc_f) && nrow(sc_f)) {
      links <- merge(cc_f, sc_f, by = c("neutral_mass_id", "formula"), suffixes = c("_broad", "_std"))

      if (nrow(links)) {
        broad_links <- stats::aggregate(
          label_std ~ broad_row,
          data = links,
          FUN = function(x) paste(sort(unique(x)), collapse = ";")
        )
        std_links <- stats::aggregate(
          label_broad ~ std_row,
          data = links,
          FUN = function(x) paste(sort(unique(x)), collapse = ";")
        )

        formula_link_broad <- broad_links$label_std[
          match(cc_left$.broad_id, broad_links$broad_row)
        ]
        formula_link_std <- std_links$label_broad[
          match(sc_left$.std_id, std_links$std_row)
        ]
      }
    }
  }

  broad_only <- broad_only_rows(cc_left)
  broad_only$possible_link_formula_only <- formula_link_broad

  standards_only <- sc_left |>
    dplyr::transmute(
      neutral_mass_id,
      source = "standards_only",
      identity_match_type = NA_character_,
      possible_link_formula_only = formula_link_std,
      candidate_ppm_error,
      candidate_neutral_mass,
      broad_db_name = NA_character_,
      broad_db_formula = NA_character_,
      broad_db_inchikey = NA_character_,
      broad_db_inchi = NA_character_,
      broad_db_smiles = NA_character_,
      broad_db_id = NA_character_,
      broad_db_source = NA_character_,
      standard_db_name = candidate_name,
      standard_db_formula = candidate_formula,
      standard_db_inchikey = candidate_inchikey,
      standard_db_smiles = candidate_smiles,
      standard_db_compound_id = candidate_db_id,
      standard_db_adducts = standard_db_adducts
    )

  dplyr::bind_rows(both_rows, broad_only, standards_only) |>
    dplyr::arrange(neutral_mass_id, candidate_ppm_error)
}


#' Standard-adduct recovery score for a single candidate row
#'
#' @description
#' Computes, for one standard-linked candidate compound, the fraction of its
#' `standards_db`-annotated adducts that were also detected computationally
#' for the same neutral mass:
#' \deqn{score = |detected \cap standard| / |standard|}
#'
#' Returns `NA` (not `0`) when there is no standard-library compound within
#' tolerance for that neutral mass, since absence of an experimental
#' reference is not negative evidence.
#'
#' When several standard-library compounds (isomers) fall within tolerance of
#' the same neutral mass, each is scored independently against the same
#' `detected_adducts`, which can help discriminate between isomers based on
#' which one's expected adducts were actually observed.
#'
#' @param standard_db_adducts Character scalar: semicolon-separated adducts
#'   annotated in `standards_db` for this candidate compound, or `NA` if this
#'   candidate has no standards-library link.
#' @param detected_adducts Character vector of adducts detected
#'   computationally for the same `neutral_mass_id`.
#'
#' @return A single numeric value in the range 0 to 1, or `NA_real_`.
#' @keywords internal
compute_standard_adduct_recovery <- function(standard_db_adducts, detected_adducts) {
  if (is.na(standard_db_adducts) || !nzchar(standard_db_adducts)) {
    return(NA_real_)
  }

  std_set <- strsplit(standard_db_adducts, ";", fixed = TRUE)[[1]]
  std_set <- unique(trimws(std_set))
  std_set <- std_set[nzchar(std_set)]

  if (!length(std_set)) {
    return(NA_real_)
  }

  det_set <- unique(detected_adducts[!is.na(detected_adducts) & nzchar(detected_adducts)])

  length(intersect(std_set, det_set)) / length(std_set)
}
