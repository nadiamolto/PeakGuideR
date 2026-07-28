#' Plot an ion image
#'
#' @description
#' Plots the spatial ion image (per-pixel intensity) of a single feature in a
#' peak matrix, using `ggplot2`.
#'
#' @param pkm Peak matrix list containing at least `mass` (numeric vector),
#'   `intensity` (numeric pixel x feature matrix) and `pos` (matrix or
#'   data.frame with `x`, `y` pixel coordinate columns).
#' @param idx Optional integer. Column index of the feature to plot in
#'   `pkm$mass`/`pkm$intensity` (the same indexing convention used by
#'   `feature_summary$idx`, `morph_results$idx_M0`/`idx_cand`,
#'   `adduct_edges$idx_i`/`idx_j`, `candidate_annotations$feature_idx` and
#'   `adduct_fam$family_members$idx`). This is the preferred way to select a
#'   feature. Exactly one of `idx` or `mz` must be supplied.
#' @param mz Optional numeric. Target m/z used to look up the closest feature
#'   in `pkm$mass` within `tol_ppm`, for exploratory use when a feature index
#'   is not already known. Exactly one of `idx` or `mz` must be supplied.
#' @param tol_ppm PPM tolerance used when searching by `mz`. Ignored (without
#'   a warning) when `idx` is supplied.
#' @param palette Colour palette for the intensity scale. `"viridis"`
#'   (default) uses the blue-to-yellow `ggplot2::scale_fill_viridis_c()`
#'   scale. Any other value is passed through as the `option` argument of
#'   `scale_fill_viridis_c()` (e.g. `"magma"`, `"inferno"`).
#' @param title Optional plot title. If `NULL` (default), a title is
#'   generated from the plotted feature's m/z - when the feature was located
#'   via `mz`, the title reflects the m/z actually found, not the one
#'   requested.
#' @param flip_y Logical. If `TRUE` (default), the pixel y coordinate is
#'   negated before plotting, matching the top-left image origin convention
#'   used elsewhere in the package.
#' @param clip_quantile Optional numeric in `(0, 1]`. If supplied, intensity
#'   values above this quantile are clipped before mapping to colour, so a
#'   handful of saturated pixels do not flatten the rest of the colour scale.
#' @param show_legend Logical. If `FALSE`, the intensity colour legend is
#'   omitted.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' data(example_pkm, package = "PeakGuideR")
#'
#' plot_ion_image(example_pkm, idx = 1)
#' plot_ion_image(example_pkm, mz = 762.6, tol_ppm = 50)
#' }
#'
#' @export
plot_ion_image <- function(
    pkm,
    idx = NULL,
    mz = NULL,
    tol_ppm = 10,
    palette = "viridis",
    title = NULL,
    flip_y = TRUE,
    clip_quantile = NULL,
    show_legend = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required to use `plot_ion_image()`. ",
      "Please install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }

  if (!is.list(pkm) || is.null(pkm$mass) || is.null(pkm$intensity) || is.null(pkm$pos)) {
    stop("`pkm` must be a list containing `mass`, `intensity` and `pos`.", call. = FALSE)
  }
  if (!is.numeric(pkm$mass)) {
    stop("`pkm$mass` must be a numeric vector.", call. = FALSE)
  }
  if (!is.matrix(pkm$intensity) || !is.numeric(pkm$intensity)) {
    stop("`pkm$intensity` must be a numeric matrix.", call. = FALSE)
  }
  if (length(pkm$mass) != ncol(pkm$intensity)) {
    stop("Length of `pkm$mass` must match the number of columns in `pkm$intensity`.", call. = FALSE)
  }
  if (!all(c("x", "y") %in% colnames(pkm$pos))) {
    stop("`pkm$pos` must have `x` and `y` columns.", call. = FALSE)
  }
  if (nrow(pkm$pos) != nrow(pkm$intensity)) {
    stop("Number of rows in `pkm$pos` must match the number of rows in `pkm$intensity`.", call. = FALSE)
  }

  if (is.null(idx) && is.null(mz)) {
    stop("Exactly one of `idx` or `mz` must be supplied.", call. = FALSE)
  }
  if (!is.null(idx) && !is.null(mz)) {
    stop("Only one of `idx` or `mz` may be supplied, not both.", call. = FALSE)
  }

  if (!is.null(idx)) {
    idx <- as.integer(idx)
    if (length(idx) != 1 || is.na(idx) || idx < 1 || idx > length(pkm$mass)) {
      stop("`idx` must be a single valid column index into `pkm$mass`/`pkm$intensity`.", call. = FALSE)
    }
  } else {
    if (!is.numeric(mz) || length(mz) != 1 || is.na(mz)) {
      stop("`mz` must be a single numeric value.", call. = FALSE)
    }
    ppm_diff <- 1e6 * abs(pkm$mass - mz) / mz
    idx <- which.min(ppm_diff)
    if (ppm_diff[idx] > tol_ppm) {
      stop(
        sprintf(
          "No feature found within %.1f ppm of m/z %.5f. Closest match: m/z %.5f (%.1f ppm away).",
          tol_ppm, mz, pkm$mass[idx], ppm_diff[idx]
        ),
        call. = FALSE
      )
    }
  }

  intensity_vals <- as.numeric(pkm$intensity[, idx])

  if (!is.null(clip_quantile)) {
    if (!is.numeric(clip_quantile) || length(clip_quantile) != 1 ||
        clip_quantile <= 0 || clip_quantile > 1) {
      stop("`clip_quantile` must be a single numeric value in (0, 1].", call. = FALSE)
    }
    clip_val <- stats::quantile(intensity_vals, probs = clip_quantile, na.rm = TRUE)
    intensity_vals <- pmin(intensity_vals, clip_val)
  }

  df <- data.frame(
    x = as.numeric(pkm$pos[, "x"]),
    y = if (isTRUE(flip_y)) -as.numeric(pkm$pos[, "y"]) else as.numeric(pkm$pos[, "y"]),
    intensity = intensity_vals
  )

  if (is.null(title)) {
    title <- sprintf("m/z %.4f", pkm$mass[idx])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$intensity)) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, fill = "Intensity") +
    ggplot2::theme_void()

  p <- if (identical(palette, "viridis")) {
    p + ggplot2::scale_fill_viridis_c(option = "D")
  } else {
    p + ggplot2::scale_fill_viridis_c(option = palette)
  }

  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}
