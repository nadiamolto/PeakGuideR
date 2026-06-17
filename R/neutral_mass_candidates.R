#' Build neutral-mass candidate table
#'
#' @description
#' Builds a neutral-mass candidate table from adduct families.
#'
#' This function summarizes neutral masses inferred from adduct families,
#' adds isotope and EIPS support from `feature_summary`, matches the inferred
#' neutral masses against a compound mass database, and optionally checks
#' whether each candidate compound and inferred adduct are supported by the
#' standard-adduct library.
#'
#' The output contains one row per inferred neutral mass and compound candidate.
#' Candidate compounds are putative mass matches, not definitive identifications.
#'
#' @param adduct_fam Output from `adduct_families()`.
#' @param feature_summary Output from `build_feature_summary()`.
#' @param compound_db Compound mass database. If `NULL`, the included
#'   non-commercial compound mass database is loaded.
#' @param standards_db Standard adduct library. If `NULL` and `matrix = "HCCA"`,
#'   the included non-commercial standard adduct library is loaded.
#' @param ion_mode `"pos"` or `"neg"`.
#' @param matrix Matrix name. Standard-adduct support is currently applied only
#'   when `matrix = "HCCA"`. Use `NULL` to skip standard-adduct support.
#' @param ppm_tol PPM tolerance for compound mass matching.
#' @param neutral_cluster_ppm PPM tolerance used to group similar neutral masses
#'   inferred from different adduct families.
#' @param top_n Maximum number of compound candidates kept per neutral mass.
#' @param quiet Logical. If `FALSE`, database loading functions may print notices.
#'
#' @return A data.frame with one row per neutral mass and compound candidate.
#' @export
build_neutral_mass_candidates <- function(
    adduct_fam,
    feature_summary,
    compound_db = NULL,
    standards_db = NULL,
    ion_mode = c("pos", "neg"),
    matrix = NULL,
    ppm_tol = 5,
    neutral_cluster_ppm = 5,
    top_n = 10L,
    quiet = FALSE
) {
  ion_mode <- match.arg(ion_mode)

  stopifnot(is.list(adduct_fam))
  stopifnot(is.data.frame(feature_summary))

  if (!all(c("family_summary", "family_members") %in% names(adduct_fam))) {
    stop(
      "adduct_fam must contain `family_summary` and `family_members`.",
      call. = FALSE
    )
  }

  stopifnot(is.data.frame(adduct_fam$family_summary))
  stopifnot(is.data.frame(adduct_fam$family_members))

  required_family_cols <- c("family_id", "neutral_mass_consensus")
  missing_family_cols <- setdiff(
    required_family_cols,
    names(adduct_fam$family_summary)
  )

  if (length(missing_family_cols) > 0) {
    stop(
      "adduct_fam$family_summary is missing required columns: ",
      paste(missing_family_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_member_cols <- c("family_id", "idx", "mz", "adduct")
  missing_member_cols <- setdiff(
    required_member_cols,
    names(adduct_fam$family_members)
  )

  if (length(missing_member_cols) > 0) {
    stop(
      "adduct_fam$family_members is missing required columns: ",
      paste(missing_member_cols, collapse = ", "),
      call. = FALSE
    )
  }

  top_n <- as.integer(top_n)
  if (!is.finite(top_n) || top_n < 1L) {
    top_n <- 10L
  }

  ppm_tol <- as.numeric(ppm_tol)
  if (!is.finite(ppm_tol) || ppm_tol <= 0) {
    ppm_tol <- 5
  }

  neutral_cluster_ppm <- as.numeric(neutral_cluster_ppm)
  if (!is.finite(neutral_cluster_ppm) || neutral_cluster_ppm <= 0) {
    neutral_cluster_ppm <- ppm_tol
  }

  if (is.null(compound_db)) {
    compound_db <- load_compound_mass_database(quiet = quiet)
  }

  stopifnot(is.data.frame(compound_db))

  if (!"MonoisotopicMass" %in% names(compound_db)) {
    stop(
      "compound_db must contain a `MonoisotopicMass` column.",
      call. = FALSE
    )
  }

  fam_summary <- adduct_fam$family_summary
  fam_members <- adduct_fam$family_members

  fam_summary$.neutral_mass <- as.numeric(fam_summary$neutral_mass_consensus)

  fam_summary <- fam_summary |>
    dplyr::filter(is.finite(.neutral_mass)) |>
    dplyr::arrange(.neutral_mass)

  if (nrow(fam_summary) == 0) {
    return(dplyr::tibble())
  }

  neutral_mass_id <- integer(nrow(fam_summary))
  current_id <- 1L
  neutral_mass_id[1] <- current_id
  current_mass <- fam_summary$.neutral_mass[1]

  if (nrow(fam_summary) > 1) {
    for (i in 2:nrow(fam_summary)) {
      ppm_diff <- 1e6 * abs(fam_summary$.neutral_mass[i] - current_mass) /
        pmax(abs(current_mass), 1e-12)

      if (is.finite(ppm_diff) && ppm_diff <= neutral_cluster_ppm) {
        neutral_mass_id[i] <- current_id
      } else {
        current_id <- current_id + 1L
        neutral_mass_id[i] <- current_id
        current_mass <- fam_summary$.neutral_mass[i]
      }
    }
  }

  fam_summary$neutral_mass_id <- neutral_mass_id

  member_info <- fam_members |>
    dplyr::left_join(
      fam_summary |>
        dplyr::select(family_id, neutral_mass_id),
      by = "family_id"
    ) |>
    dplyr::filter(!is.na(neutral_mass_id))

  feature_cols <- intersect(
    c(
      "idx",
      "is_c13_m0",
      "c13_score",
      "has_c13_m2_support",
      "c13_m2_score",
      "has_eips",
      "eips_elements",
      "eips_score"
    ),
    names(feature_summary)
  )

  member_feature_info <- member_info |>
    dplyr::left_join(
      feature_summary |>
        dplyr::select(dplyr::all_of(feature_cols)),
      by = "idx"
    ) |>
    dplyr::mutate(
      is_c13_m0 = dplyr::coalesce(is_c13_m0, FALSE),
      has_c13_m2_support = dplyr::coalesce(has_c13_m2_support, FALSE),
      has_eips = dplyr::coalesce(has_eips, FALSE)
    )

  neutral_summary <- fam_summary |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      neutral_mass_consensus = mean(.neutral_mass, na.rm = TRUE),
      neutral_mass_min = min(.neutral_mass, na.rm = TRUE),
      neutral_mass_max = max(.neutral_mass, na.rm = TRUE),
      neutral_mass_range_ppm =
        1e6 * (neutral_mass_max - neutral_mass_min) /
        pmax(abs(neutral_mass_consensus), 1e-12),
      n_adduct_families = dplyr::n_distinct(family_id),
      adduct_family_ids = paste(sort(unique(family_id)), collapse = ";"),
      .groups = "drop"
    )

  adduct_support_summary <- member_feature_info |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::summarise(
      n_features = dplyr::n_distinct(idx),
      feature_idx = paste(sort(unique(idx)), collapse = ";"),
      feature_mz = paste(round(sort(unique(mz)), 6), collapse = ";"),
      n_adducts_inferred = dplyr::n_distinct(adduct),
      inferred_adducts = paste(sort(unique(adduct)), collapse = ";"),

      has_CIR_support = any(is_c13_m0, na.rm = TRUE),
      CIR_feature_idx = paste(
        sort(unique(idx[is_c13_m0 %in% TRUE])),
        collapse = ";"
      ),
      CIR_score = if (any(is.finite(c13_score))) {
        max(c13_score, na.rm = TRUE)
      } else {
        NA_real_
      },

      has_C13_M2_support = any(has_c13_m2_support, na.rm = TRUE),
      C13_M2_score = if (any(is.finite(c13_m2_score))) {
        max(c13_m2_score, na.rm = TRUE)
      } else {
        NA_real_
      },

      has_EIPS_support = any(has_eips, na.rm = TRUE),
      EIPS_elements = paste(
        sort(unique(unlist(strsplit(
          paste(eips_elements[has_eips %in% TRUE], collapse = ";"),
          ";",
          fixed = TRUE
        )))),
        collapse = ";"
      ),
      EIPS_score = if (any(is.finite(eips_score))) {
        max(eips_score, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      CIR_feature_idx = dplyr::if_else(
        CIR_feature_idx == "",
        NA_character_,
        CIR_feature_idx
      ),
      EIPS_elements = dplyr::if_else(
        EIPS_elements == "",
        NA_character_,
        EIPS_elements
      )
    )

  neutral_table <- neutral_summary |>
    dplyr::left_join(adduct_support_summary, by = "neutral_mass_id")

  db_cols <- intersect(
    c(
      "Source",
      "DB_ID",
      "Name",
      "MolecularFormula",
      "MonoisotopicMass",
      "StdInChI",
      "StdInChIKey",
      "SMILES",
      "Kegg"
    ),
    names(compound_db)
  )

  db2 <- compound_db |>
    dplyr::select(dplyr::all_of(db_cols)) |>
    dplyr::mutate(
      candidate_mass = as.numeric(MonoisotopicMass)
    ) |>
    dplyr::filter(is.finite(candidate_mass))

  neutral_for_matching <- neutral_table |>
    dplyr::select(neutral_mass_id, neutral_mass_consensus)

  compound_candidates <- tidyr::crossing(neutral_for_matching, db2) |>
    dplyr::mutate(
      candidate_ppm_error =
        1e6 * abs(candidate_mass - neutral_mass_consensus) /
        pmax(abs(neutral_mass_consensus), 1e-12)
    ) |>
    dplyr::filter(
      is.finite(candidate_ppm_error),
      candidate_ppm_error <= ppm_tol
    ) |>
    dplyr::arrange(neutral_mass_id, candidate_ppm_error) |>
    dplyr::group_by(neutral_mass_id) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup() |>
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

  if (nrow(compound_candidates) == 0) {
    out <- neutral_table |>
      dplyr::mutate(
        candidate_source = NA_character_,
        candidate_db_id = NA_character_,
        candidate_name = NA_character_,
        candidate_formula = NA_character_,
        candidate_neutral_mass = NA_real_,
        candidate_inchi = NA_character_,
        candidate_inchikey = NA_character_,
        candidate_smiles = NA_character_,
        candidate_kegg = NA_character_,
        candidate_mass = NA_real_,
        candidate_ppm_error = NA_real_,
      )
  } else {
    out <- neutral_table |>
      dplyr::left_join(
        compound_candidates,
        by = c("neutral_mass_id", "neutral_mass_consensus")
      )
  }

  out <- out |>
    dplyr::mutate(
      matrix = if (is.null(matrix)) NA_character_ else as.character(matrix),
      has_standard_compound_match = NA,
      has_standard_adduct_match = NA,
      standard_matched_name = NA_character_,
      standard_matched_adducts = NA_character_,
      standard_matched_inferred_adducts = NA_character_,
      n_standard_matched_adducts = NA_integer_
    )

  if (!is.null(matrix) && identical(toupper(matrix), "HCCA")) {
    if (is.null(standards_db)) {
      standards_db <- load_standards_adduct_library(quiet = quiet)
    }

    stopifnot(is.data.frame(standards_db))

    required_std_cols <- c("adduct", "POLARITY")
    missing_std_cols <- setdiff(required_std_cols, names(standards_db))

    if (length(missing_std_cols) > 0) {
      warning(
        "standards_db is missing required columns: ",
        paste(missing_std_cols, collapse = ", "),
        ". Standard-adduct support will not be added.",
        call. = FALSE
      )
    } else {
      std <- standards_db |>
        dplyr::filter(POLARITY == ion_mode)

      if ("found" %in% names(std)) {
        std <- std |>
          dplyr::filter(found %in% TRUE)
      }

      if ("keep_final" %in% names(std)) {
        std <- std |>
          dplyr::filter(keep_final %in% TRUE)
      }

      if (!"COMPOUND_ID" %in% names(std)) {
        std$COMPOUND_ID <- NA_character_
      }

      if (!"Master_List_NAME" %in% names(std)) {
        std$Master_List_NAME <- NA_character_
      }

      if (!"HMDB_clean" %in% names(std)) {
        std$HMDB_clean <- NA_character_
      }

      if (!"ChEBI" %in% names(std)) {
        std$ChEBI <- NA_character_
      }

      if (!"InCHIKey" %in% names(std)) {
        std$InCHIKey <- NA_character_
      }

      if (!"SMILES" %in% names(std)) {
        std$SMILES <- NA_character_
      }

      if (!"MOLECULAR_FORMULA" %in% names(std)) {
        std$MOLECULAR_FORMULA <- NA_character_
      }

      std_clean <- std |>
        dplyr::mutate(
          std_compound_id = as.character(COMPOUND_ID),
          std_name = as.character(Master_List_NAME),
          std_hmdb = as.character(HMDB_clean),
          std_chebi = as.character(ChEBI),
          std_inchikey = as.character(InCHIKey),
          std_smiles = as.character(SMILES),
          std_formula = as.character(MOLECULAR_FORMULA),
          std_adduct = as.character(adduct)
        ) |>
        dplyr::select(
          std_compound_id,
          std_name,
          std_hmdb,
          std_chebi,
          std_inchikey,
          std_smiles,
          std_formula,
          std_adduct
        ) |>
        dplyr::distinct()

      candidate_adducts <- out |>
        dplyr::select(
          neutral_mass_id,
          candidate_source,
          candidate_db_id,
          candidate_name,
          candidate_formula,
          candidate_inchikey,
          candidate_smiles,
          inferred_adducts
        ) |>
        dplyr::mutate(
          inferred_adducts = as.character(inferred_adducts)
        ) |>
        tidyr::separate_rows(inferred_adducts, sep = ";") |>
        dplyr::rename(inferred_adduct = inferred_adducts) |>
        dplyr::filter(
          !is.na(candidate_name),
          !is.na(inferred_adduct),
          inferred_adduct != ""
        )

      candidate_adducts <- candidate_adducts |>
        dplyr::mutate(
          match_hmdb = dplyr::if_else(
            candidate_source == "HMDB",
            as.character(candidate_db_id),
            NA_character_
          ),
          match_chebi = dplyr::if_else(
            candidate_source == "CHEBI",
            as.character(candidate_db_id),
            NA_character_
          ),
          match_inchikey = as.character(candidate_inchikey),
          match_smiles = as.character(candidate_smiles),
          match_formula = as.character(candidate_formula)
        )

      std_matches <- dplyr::bind_rows(
        candidate_adducts |>
          dplyr::filter(!is.na(match_inchikey), match_inchikey != "") |>
          dplyr::inner_join(
            std_clean |>
              dplyr::filter(!is.na(std_inchikey), std_inchikey != ""),
            by = c("match_inchikey" = "std_inchikey"),
            relationship = "many-to-many"
          ),

        candidate_adducts |>
          dplyr::filter(!is.na(match_hmdb), match_hmdb != "") |>
          dplyr::inner_join(
            std_clean |>
              dplyr::filter(!is.na(std_hmdb), std_hmdb != ""),
            by = c("match_hmdb" = "std_hmdb"),
            relationship = "many-to-many"
          ),

        candidate_adducts |>
          dplyr::filter(!is.na(match_chebi), match_chebi != "") |>
          dplyr::inner_join(
            std_clean |>
              dplyr::filter(!is.na(std_chebi), std_chebi != ""),
            by = c("match_chebi" = "std_chebi"),
            relationship = "many-to-many"
          ),

        candidate_adducts |>
          dplyr::filter(!is.na(match_smiles), match_smiles != "") |>
          dplyr::inner_join(
            std_clean |>
              dplyr::filter(!is.na(std_smiles), std_smiles != ""),
            by = c("match_smiles" = "std_smiles"),
            relationship = "many-to-many"
          ),

        candidate_adducts |>
          dplyr::filter(!is.na(match_formula), match_formula != "") |>
          dplyr::inner_join(
            std_clean |>
              dplyr::filter(!is.na(std_formula), std_formula != ""),
            by = c("match_formula" = "std_formula"),
            relationship = "many-to-many"
          )
      ) |>
        dplyr::distinct()

      if (nrow(std_matches) > 0) {
        std_support <- std_matches |>
          dplyr::mutate(
            inferred_adduct_is_in_standard =
              inferred_adduct == std_adduct
          ) |>
          dplyr::group_by(
            neutral_mass_id,
            candidate_source,
            candidate_db_id,
            candidate_name,
            candidate_formula
          ) |>
          dplyr::summarise(
            has_standard_compound_match = TRUE,
            has_standard_adduct_match = any(
              inferred_adduct_is_in_standard,
              na.rm = TRUE
            ),
            standard_matched_name = paste(
              sort(unique(std_name)),
              collapse = ";"
            ),
            standard_matched_adducts = paste(
              sort(unique(std_adduct)),
              collapse = ";"
            ),
            standard_matched_inferred_adducts = paste(
              sort(unique(inferred_adduct[inferred_adduct_is_in_standard])),
              collapse = ";"
            ),
            n_standard_matched_adducts = dplyr::n_distinct(
              inferred_adduct[inferred_adduct_is_in_standard]
            ),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            standard_matched_inferred_adducts = dplyr::if_else(
              standard_matched_inferred_adducts == "",
              NA_character_,
              standard_matched_inferred_adducts
            )
          )

        out <- out |>
          dplyr::select(
            -has_standard_compound_match,
            -has_standard_adduct_match,
            -standard_matched_name,
            -standard_matched_adducts,
            -standard_matched_inferred_adducts,
            -n_standard_matched_adducts
          ) |>
          dplyr::left_join(
            std_support,
            by = c(
              "neutral_mass_id",
              "candidate_source",
              "candidate_db_id",
              "candidate_name",
              "candidate_formula"
            )
          ) |>
          dplyr::mutate(
            has_standard_compound_match = dplyr::coalesce(
              has_standard_compound_match,
              FALSE
            ),
            has_standard_adduct_match = dplyr::coalesce(
              has_standard_adduct_match,
              FALSE
            ),
            n_standard_matched_adducts = dplyr::coalesce(
              as.integer(n_standard_matched_adducts),
              0L
            )
          )
      } else {
        out <- out |>
          dplyr::mutate(
            has_standard_compound_match = FALSE,
            has_standard_adduct_match = FALSE,
            n_standard_matched_adducts = 0L
          )
      }
    }
  }

  out <- out |>
    dplyr::select(
      -dplyr::any_of(c(
        "neutral_mass_min",
        "neutral_mass_max",
        "neutral_mass_range_ppm"
      ))
    )

  out |>
    dplyr::arrange(
      neutral_mass_id,
      candidate_ppm_error
    )
}
