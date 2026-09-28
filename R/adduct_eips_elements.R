#' EIPS-tracked element introduced by an adduct
#'
#' @description
#' Maps an adduct name (as in `default_adducts()$name`) to the elemental
#' isotope-pattern (EIPS) element it adds to the ion, if any. EIPS tracks
#' N, O, S, Cl and Br; only adducts that carry one of these into the ion have
#' an entry.
#'
#' In the current [default_adducts()] table `[M+Cl]-` is the only such adduct
#' (H, Na, K and NH4 add none of the tracked elements), so a `35Cl`/`37Cl`
#' isotope pattern on a feature annotated as `[M+Cl]-` is explained by the
#' adduct itself and is not evidence that the neutral molecule contains
#' chlorine. This is a fixed internal lookup for the default adduct table, not
#' a formula parser.
#'
#' @param adduct Character vector of adduct names.
#'
#' @return Character vector the same length as `adduct`: the introduced element,
#'   or `NA_character_` where the adduct introduces none.
#'
#' @keywords internal
adduct_introduced_eips_element <- function(adduct) {
  lookup <- c("[M+Cl]-" = "Cl")
  unname(lookup[as.character(adduct)])
}


#' Drop adduct-explained elements from per-feature EIPS evidence
#'
#' @description
#' `iso_morphology_candidates()` / [eips_score()] detect elemental isotope
#' patterns from mass difference and spatial correlation alone, with no
#' knowledge of the adduct assigned to the same feature elsewhere in the
#' pipeline. When the detected element is one the assigned adduct itself
#' introduces (see [adduct_introduced_eips_element()]), the pattern comes from
#' the adduct, not the neutral molecule, and must not count as EIPS evidence
#' for the neutral-mass candidate.
#'
#' This removes such elements from `eips_elements` on a per-feature basis. A
#' feature left with no elements loses its EIPS support (`has_eips = FALSE`,
#' `eips_score = NA`). Evidence for elements the adduct does not introduce -
#' any N/O/S evidence, or Cl/Br evidence on a feature assigned a different
#' adduct - is returned unchanged.
#'
#' @param adduct Character vector of assigned adduct names, one per feature.
#' @param has_eips Logical vector, the feature's raw EIPS flag.
#' @param eips_elements Character vector of `";"`-separated EIPS elements, `NA`
#'   where the feature has no EIPS evidence.
#' @param eips_score Numeric vector, the feature's raw EIPS score.
#'
#' @return A data.frame with adjusted `has_eips`, `eips_elements` and
#'   `eips_score`, one row per input feature.
#'
#' @keywords internal
drop_adduct_explained_eips <- function(adduct, has_eips, eips_elements, eips_score) {
  has_eips <- dplyr::coalesce(as.logical(has_eips), FALSE)
  eips_elements <- as.character(eips_elements)
  eips_score <- as.numeric(eips_score)
  introduced <- adduct_introduced_eips_element(adduct)

  kept_elements <- character(length(has_eips))
  still_has_eips <- logical(length(has_eips))

  for (i in seq_along(has_eips)) {
    if (!isTRUE(has_eips[i])) {
      kept_elements[i] <- NA_character_
      still_has_eips[i] <- FALSE
      next
    }
    # EIPS present but the detected element is unknown, or the adduct adds no
    # tracked element: nothing can be attributed to the adduct, keep as-is.
    if (is.na(eips_elements[i]) || is.na(introduced[i])) {
      kept_elements[i] <- eips_elements[i]
      still_has_eips[i] <- TRUE
      next
    }
    els <- trimws(strsplit(eips_elements[i], ";", fixed = TRUE)[[1]])
    els <- setdiff(els[nzchar(els)], introduced[i])
    if (length(els) == 0) {
      kept_elements[i] <- NA_character_
      still_has_eips[i] <- FALSE
    } else {
      kept_elements[i] <- paste(els, collapse = ";")
      still_has_eips[i] <- TRUE
    }
  }

  # eips_score is a per-feature max over EIPS relations and cannot be split by
  # element here; it is kept as-is while any element survives (so genuine
  # evidence scores exactly as before) and cleared only when none does.
  data.frame(
    has_eips = still_has_eips,
    eips_elements = kept_elements,
    eips_score = dplyr::if_else(still_has_eips, eips_score, NA_real_),
    stringsAsFactors = FALSE
  )
}
