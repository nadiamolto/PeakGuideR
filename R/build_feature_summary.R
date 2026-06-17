#' Build feature-level summary from PeakGuideR relations
#'
#' @description
#' Builds a feature-level summary from a PeakGuideR relation table.
#'
#' The output contains one row per feature and summarizes whether each peak is
#' involved in isotope, EIPS or adduct relationships.
#'
#' @param relation_table Output from `build_relation_table()`.
#' @param adduct_fam Optional output from `adduct_families()`.
#' @param pkm Optional peak matrix object with `mass`. If supplied, all features
#'   in `pkm$mass` are included. If `NULL`, only features present in
#'   `relation_table` are returned.
#'
#' @return A data.frame with one row per feature.
#' @export
build_feature_summary <- function(
    relation_table,
    adduct_fam = NULL,
    pkm = NULL
) {
  stopifnot(is.data.frame(relation_table))

  required_cols <- c(
    "from_idx", "to_idx",
    "from_mz", "to_mz",
    "relation_type", "evidence_type",
    "evidence_score", "from_role", "to_role"
  )

  missing_cols <- setdiff(required_cols, names(relation_table))
  if (length(missing_cols) > 0) {
    stop(
      "relation_table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(pkm)) {
    stopifnot(is.list(pkm), !is.null(pkm$mass))

    features <- data.frame(
      idx = seq_along(pkm$mass),
      mz = as.numeric(pkm$mass),
      stringsAsFactors = FALSE
    )
  } else {
    features <- dplyr::bind_rows(
      relation_table |>
        dplyr::transmute(idx = from_idx, mz = from_mz),
      relation_table |>
        dplyr::transmute(idx = to_idx, mz = to_mz)
    ) |>
      dplyr::filter(!is.na(idx), !is.na(mz)) |>
      dplyr::distinct(idx, .keep_all = TRUE)
  }

  c13_m0 <- relation_table |>
    dplyr::filter(
      evidence_type == "CIR",
      relation_type == "C13_M1"
    ) |>
    dplyr::group_by(idx = from_idx) |>
    dplyr::slice_max(
      order_by = evidence_score,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx,
      is_c13_m0 = TRUE,
      c13_m1_idx = to_idx,
      c13_m1_mz = to_mz,
      c13_score = evidence_score
    )

  c13_m1 <- relation_table |>
    dplyr::filter(
      evidence_type == "CIR",
      relation_type == "C13_M1"
    ) |>
    dplyr::group_by(idx = to_idx) |>
    dplyr::slice_max(
      order_by = evidence_score,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx,
      is_c13_m1 = TRUE,
      c13_parent_idx = from_idx,
      c13_parent_mz = from_mz,
      c13_parent_score = evidence_score
    )

  c13_m2_support <- relation_table |>
    dplyr::filter(
      evidence_type == "isotope_morphology",
      relation_type == "C13_M2"
    ) |>
    dplyr::group_by(idx = from_idx) |>
    dplyr::slice_max(
      order_by = evidence_score,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx,
      has_c13_m2_support = TRUE,
      c13_m2_idx = to_idx,
      c13_m2_mz = to_mz,
      c13_m2_score = evidence_score
    )

  c13_m2 <- relation_table |>
    dplyr::filter(
      evidence_type == "isotope_morphology",
      relation_type == "C13_M2"
    ) |>
    dplyr::group_by(idx = to_idx) |>
    dplyr::slice_max(
      order_by = evidence_score,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx,
      is_c13_m2 = TRUE,
      c13_m2_parent_idx = from_idx,
      c13_m2_parent_mz = from_mz,
      c13_m2_parent_score = evidence_score
    )

  eips_parent <- relation_table |>
    dplyr::filter(grepl("^EIPS", evidence_type)) |>
    dplyr::group_by(idx = from_idx) |>
    dplyr::summarise(
      has_eips = TRUE,
      eips_elements = paste(
        sort(unique(gsub("^EIPS_", "", evidence_type))),
        collapse = ";"
      ),
      eips_partner_idx = paste(sort(unique(to_idx)), collapse = ";"),
      eips_score = max(evidence_score, na.rm = TRUE),
      .groups = "drop"
    )

  eips_isotope <- relation_table |>
    dplyr::filter(grepl("^EIPS", evidence_type)) |>
    dplyr::group_by(idx = to_idx) |>
    dplyr::slice_max(
      order_by = evidence_score,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx,
      is_eips_isotope = TRUE,
      eips_parent_idx = from_idx,
      eips_parent_mz = from_mz,
      eips_relation_type = relation_type,
      eips_parent_score = evidence_score
    )

  adduct_roles_long <- dplyr::bind_rows(
    relation_table |>
      dplyr::filter(evidence_type == "adduct_family") |>
      dplyr::transmute(
        idx = from_idx,
        adduct_family_id = group_id,
        adduct_role = from_role,
        adduct_score = evidence_score
      ),
    relation_table |>
      dplyr::filter(evidence_type == "adduct_family") |>
      dplyr::transmute(
        idx = to_idx,
        adduct_family_id = group_id,
        adduct_role = to_role,
        adduct_score = evidence_score
      )
  ) |>
    dplyr::filter(!is.na(idx))

  adduct_summary <- adduct_roles_long |>
    dplyr::group_by(idx) |>
    dplyr::summarise(
      has_adduct_family = TRUE,
      adduct_family_ids = paste(
        sort(unique(adduct_family_id)),
        collapse = ";"
      ),
      n_adduct_families = dplyr::n_distinct(adduct_family_id),
      adduct_roles = paste(
        sort(unique(adduct_role)),
        collapse = ";"
      ),
      n_adduct_roles = dplyr::n_distinct(adduct_role),
      adduct_best_score = max(adduct_score, na.rm = TRUE),
      adduct_role_status = dplyr::case_when(
        dplyr::n_distinct(adduct_role) == 1 ~ "single",
        dplyr::n_distinct(adduct_role) > 1 ~ "multiple",
        TRUE ~ "none"
      ),
      adduct_family_status = dplyr::case_when(
        dplyr::n_distinct(adduct_family_id) == 1 ~ "single_family",
        dplyr::n_distinct(adduct_family_id) > 1 ~ "multiple_families",
        TRUE ~ "none"
      ),
      .groups = "drop"
    )

  adduct_detail <- NULL

  if (!is.null(adduct_fam) &&
      is.list(adduct_fam) &&
      "family_members" %in% names(adduct_fam) &&
      "family_summary" %in% names(adduct_fam) &&
      is.data.frame(adduct_fam$family_members) &&
      is.data.frame(adduct_fam$family_summary)) {

    adduct_detail <- adduct_fam$family_members |>
      dplyr::left_join(
        adduct_fam$family_summary |>
          dplyr::select(
            family_id,
            neutral_mass_consensus,
            family_size,
            n_edges,
            mean_score_adduct,
            has_role_conflict
          ) |>
          dplyr::rename(
            adduct_family_neutral_mass_consensus = neutral_mass_consensus,
            adduct_family_size = family_size,
            adduct_family_n_edges = n_edges,
            adduct_family_mean_score = mean_score_adduct,
            adduct_family_has_role_conflict = has_role_conflict
          ),
        by = "family_id"
      ) |>
      dplyr::group_by(idx) |>
      dplyr::slice_max(
        order_by = mean_edge_score,
        n = 1,
        with_ties = FALSE
      ) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        idx,
        main_adduct_family_id = family_id,
        main_adduct_role = adduct,
        neutral_mass_consensus = adduct_family_neutral_mass_consensus,
        adduct_family_size,
        adduct_family_n_edges,
        adduct_family_mean_score,
        adduct_family_has_role_conflict
      )
  }

  out <- features |>
    dplyr::left_join(c13_m0, by = "idx") |>
    dplyr::left_join(c13_m1, by = "idx") |>
    dplyr::left_join(c13_m2_support, by = "idx") |>
    dplyr::left_join(c13_m2, by = "idx") |>
    dplyr::left_join(eips_parent, by = "idx") |>
    dplyr::left_join(eips_isotope, by = "idx") |>
    dplyr::left_join(adduct_summary, by = "idx")

  if (!is.null(adduct_detail)) {
    out <- out |>
      dplyr::left_join(adduct_detail, by = "idx")
  }

  out |>
    dplyr::mutate(
      is_c13_m0 = dplyr::coalesce(is_c13_m0, FALSE),
      is_c13_m1 = dplyr::coalesce(is_c13_m1, FALSE),
      has_c13_m2_support = dplyr::coalesce(has_c13_m2_support, FALSE),
      is_c13_m2 = dplyr::coalesce(is_c13_m2, FALSE),
      has_eips = dplyr::coalesce(has_eips, FALSE),
      is_eips_isotope = dplyr::coalesce(is_eips_isotope, FALSE),
      has_adduct_family = dplyr::coalesce(has_adduct_family, FALSE),

      adduct_role_status = dplyr::coalesce(adduct_role_status, "none"),
      adduct_family_status = dplyr::coalesce(adduct_family_status, "none"),

      n_adduct_families = dplyr::coalesce(as.integer(n_adduct_families), 0L),
      n_adduct_roles = dplyr::coalesce(as.integer(n_adduct_roles), 0L)
    ) |>
    dplyr::arrange(
      dplyr::desc(
        is_c13_m0 |
          is_c13_m1 |
          is_c13_m2 |
          has_c13_m2_support |
          has_eips |
          is_eips_isotope |
          has_adduct_family
      ),
      idx
    )
}
