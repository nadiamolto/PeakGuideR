#' Plot all members of an adduct family
#'
#' @description
#' Plots the ion image of every feature belonging to a given adduct family
#' (as returned by `adduct_families()`), reusing `plot_ion_image()` for each
#' panel and combining them with `patchwork`.
#'
#' @param pkm Peak matrix list, see `plot_ion_image()`.
#' @param family_id Adduct family ID to plot (matches
#'   `adduct_fam$family_members$family_id` /
#'   `adduct_fam$family_summary$family_id`).
#' @param adduct_fam List returned by `adduct_families()`, containing
#'   `family_members` and `family_summary`.
#' @param ncol Optional integer. Number of columns used to arrange the
#'   panels. If `NULL` (default), a reasonable number of columns is chosen
#'   automatically from the family size (`ceiling(sqrt(n))`).
#' @param ... Additional arguments passed through to `plot_ion_image()` (e.g.
#'   `palette`, `clip_quantile`, `flip_y`, `show_legend`).
#'
#' @return A `patchwork` object combining one ion image panel per family
#'   member.
#'
#' @examples
#' \dontrun{
#' data(example_pkm, package = "PeakGuideR")
#'
#' adduct_res <- adduct_candidates(example_pkm, ion_mode = "pos")
#' adduct_fam <- adduct_families(adduct_res)
#'
#' plot_adduct_family(
#'   example_pkm,
#'   family_id = adduct_fam$family_summary$family_id[1],
#'   adduct_fam = adduct_fam
#' )
#' }
#'
#' @export
plot_adduct_family <- function(
    pkm,
    family_id,
    adduct_fam,
    ncol = NULL,
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required to use `plot_adduct_family()`. ",
      "Please install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required to use `plot_adduct_family()`. ",
      "Please install it with install.packages('patchwork').",
      call. = FALSE
    )
  }

  if (!is.list(adduct_fam) ||
      !all(c("family_members", "family_summary") %in% names(adduct_fam)) ||
      !is.data.frame(adduct_fam$family_members) ||
      !is.data.frame(adduct_fam$family_summary)) {
    stop(
      "`adduct_fam` must be the list returned by `adduct_families()`, ",
      "containing `family_members` and `family_summary` data.frames.",
      call. = FALSE
    )
  }

  required_member_cols <- c("idx", "mz", "family_id", "adduct")
  missing_member_cols <- setdiff(required_member_cols, names(adduct_fam$family_members))
  if (length(missing_member_cols) > 0) {
    stop(
      "`adduct_fam$family_members` is missing required columns: ",
      paste(missing_member_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_summary_cols <- c("family_id", "neutral_mass_consensus")
  missing_summary_cols <- setdiff(required_summary_cols, names(adduct_fam$family_summary))
  if (length(missing_summary_cols) > 0) {
    stop(
      "`adduct_fam$family_summary` is missing required columns: ",
      paste(missing_summary_cols, collapse = ", "),
      call. = FALSE
    )
  }

  target_family_id <- family_id
  members <- adduct_fam$family_members |>
    dplyr::filter(.data$family_id == target_family_id)

  if (nrow(members) == 0) {
    stop(
      sprintf("`family_id` %s not found in `adduct_fam$family_members`.", target_family_id),
      call. = FALSE
    )
  }

  summary_row <- adduct_fam$family_summary |>
    dplyr::filter(.data$family_id == target_family_id)

  panels <- vector("list", nrow(members))
  for (i in seq_len(nrow(members))) {
    panels[[i]] <- plot_ion_image(
      pkm,
      idx = members$idx[i],
      title = sprintf("%s\nm/z %.4f", members$adduct[i], members$mz[i]),
      ...
    )
  }

  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(length(panels)))
  }

  combined <- patchwork::wrap_plots(panels, ncol = ncol)

  title_txt <- if (nrow(summary_row) > 0) {
    sprintf(
      "Adduct family %s \u2014 neutral mass %.4f",
      target_family_id, summary_row$neutral_mass_consensus[1]
    )
  } else {
    sprintf("Adduct family %s", target_family_id)
  }

  combined + patchwork::plot_annotation(title = title_txt)
}
