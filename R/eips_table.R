#' EIPS formula-support lookup table
#'
#' Database-derived lookup table used by `eips_score()` to evaluate optional
#' formula support for elemental isotope-pattern evidence.
#'
#' This table links candidate formulas and monoisotopic masses to the number of
#' atoms of selected non-carbon elements. It is used to check whether an inferred
#' neutral mass and estimated element count are compatible with database-derived
#' formula candidates.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{formula}{Molecular formula.}
#'   \item{mz_mono}{Monoisotopic mass or monoisotopic m/z reference value.}
#'   \item{element}{Element considered for isotope-pattern support: Br, Cl, N, O or S.}
#'   \item{model}{Model used for the isotope-pattern calculation.}
#'   \item{k}{Expected isotope peak order or isotope-shift multiplier.}
#'   \item{delta}{Expected m/z shift from the monoisotopic peak.}
#'   \item{n_el}{Number of atoms of the selected element in the formula.}
#'   \item{R_theo}{Expected theoretical isotope ratio for the formula-level entry.}
#' }
#'
#' @details
#' This object is used by `eips_score()` only when formula support is evaluated.
#' It is joined by `element` and filtered by neutral-mass agreement and agreement
#' between the inferred element count and the formula-derived element count.
#'
#' @source Internal PeakGuideR reference table derived from the compound mass
#' database distributed with PeakGuideR.
#'
#' @keywords datasets
#' @name eips_table
#' @docType data
#' @usage data(eips_table)
NULL


#' EIPS theoretical element-count table
#'
#' Precomputed theoretical lookup table used by `eips_score()` to compare
#' observed elemental isotope-pattern ratios against expected ratios for
#' different atom counts.
#'
#' This table contains theoretical values for Br, Cl, N, O and S. It is used to
#' infer the most compatible number of atoms of a given element from the observed
#' isotope-pattern ratio.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{element}{Element considered for isotope-pattern support: Br, Cl, N, O or S.}
#'   \item{model}{Model used for the isotope-pattern calculation.}
#'   \item{k}{Expected isotope peak order or isotope-shift multiplier.}
#'   \item{delta}{Expected m/z shift from the monoisotopic peak.}
#'   \item{p}{Isotopic abundance or model probability used in the calculation.}
#'   \item{n}{Number of atoms of the selected element.}
#'   \item{R_theo}{Expected theoretical isotope ratio for the given element and atom count.}
#' }
#'
#' @details
#' In this object, the `n` in `eips_n_table` refers to the number of atoms of
#' the selected element, not specifically to nitrogen. The table contains
#' theoretical element-count reference values for Br, Cl, N, O and S.
#'
#' @source Internal PeakGuideR theoretical reference table.
#'
#' @keywords datasets
#' @name eips_n_table
#' @docType data
#' @usage data(eips_n_table)
NULL
