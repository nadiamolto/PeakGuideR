#' Build a unified peak relation table
#'
#' @description
#' Builds a unified table of peak-to-peak relationships detected by PeakGuideR.
#'
#' This table integrates internal evidence from C13 isotope validation,
#' C13 M+2 support, elemental isotope-pattern support and adduct-family
#' relationships.
#'
#' Each relation keeps the original score from the evidence layer that generated
#' it. Therefore, `evidence_score` should be interpreted together with
#' `score_type`.
#'
#' @param cir_results Output from `cir_score()`.
#' @param eips_results Output from `eips_score()`.
#' @param adduct_fam Output from `adduct_families()`.
#' @param only_valid Logical. If `TRUE`, keeps only valid relations.
#'
#' @return A data.frame with one row per peak-to-peak relation.
#' @export
build_relation_table <- function(
    cir_results = NULL,
    eips_results = NULL,
    adduct_fam = NULL,
    only_valid = TRUE
) {
  edges <- list()

  # 1. C13 M+1 relations from CIR
  if (!is.null(cir_results) &&
      is.data.frame(cir_results) &&
      nrow(cir_results) > 0) {

    cir_edges <- cir_results |>
      dplyr::transmute(
        from_idx = idx_M0,
        to_idx = idx_M1,
        from_mz = mz_M0,
        to_mz = mz_M1,
        relation_type = "C13_M1",
        evidence_type = "CIR",
        evidence_score = cir_score,
        score_type = "isotope_ratio_score",
        is_valid = is_valid_c13,
        group_id = NA_integer_,
        from_role = "M0",
        to_role = "M+1",
        label = "C13 M+1"
      )

    edges[["cir"]] <- cir_edges

    # Optional C13 M+2 support
    if (all(c("has_C13_M2", "idx_C13_M2", "mz_C13_M2", "score_C13_M2") %in%
            names(cir_results))) {

      c13_m2_edges <- cir_results |>
        dplyr::filter(has_C13_M2 %in% TRUE) |>
        dplyr::transmute(
          from_idx = idx_M0,
          to_idx = idx_C13_M2,
          from_mz = mz_M0,
          to_mz = mz_C13_M2,
          relation_type = "C13_M2",
          evidence_type = "isotope_morphology",
          evidence_score = score_C13_M2,
          score_type = "spatial_isotope_score",
          is_valid = TRUE,
          group_id = NA_integer_,
          from_role = "M0",
          to_role = "M+2",
          label = "C13 M+2"
        )

      edges[["c13_m2"]] <- c13_m2_edges
    }
  }

  # 2. EIPS relations
  if (!is.null(eips_results) &&
      is.data.frame(eips_results) &&
      nrow(eips_results) > 0) {

    eips_edges <- eips_results |>
      dplyr::transmute(
        from_idx = idx_M0,
        to_idx = idx_cand,
        from_mz = mz_M0,
        to_mz = mz_iso,
        relation_type = iso_type,
        evidence_type = paste0("EIPS_", element),
        evidence_score = score_eips,
        score_type = "elemental_ratio_score",
        is_valid = is_valid_eips,
        group_id = NA_integer_,
        from_role = "M0",
        to_role = paste0(element, "_isotope"),
        label = paste0("EIPS ", element)
      )

    edges[["eips"]] <- eips_edges
  }

  # 3. Adduct relations
  if (!is.null(adduct_fam) &&
      is.list(adduct_fam) &&
      "family_edges" %in% names(adduct_fam) &&
      is.data.frame(adduct_fam$family_edges) &&
      nrow(adduct_fam$family_edges) > 0) {

    adduct_edges <- adduct_fam$family_edges |>
      dplyr::transmute(
        from_idx = idx_i,
        to_idx = idx_j,
        from_mz = mz_i,
        to_mz = mz_j,
        relation_type = paste(adduct_i, adduct_j, sep = " / "),
        evidence_type = "adduct_family",
        evidence_score = score_adduct,
        score_type = "spatial_adduct_score",
        is_valid = is_valid_adduct,
        group_id = family_id,
        from_role = adduct_i,
        to_role = adduct_j,
        label = paste0("Adduct family ", family_id)
      )

    edges[["adducts"]] <- adduct_edges
  }

  if (!length(edges)) {
    return(dplyr::tibble())
  }

  out <- dplyr::bind_rows(edges) |>
    dplyr::filter(!is.na(from_idx), !is.na(to_idx))

  if (isTRUE(only_valid) && "is_valid" %in% names(out)) {
    out <- out |>
      dplyr::filter(is_valid %in% TRUE)
  }

  out |>
    dplyr::arrange(from_idx, evidence_type, dplyr::desc(evidence_score))
}
