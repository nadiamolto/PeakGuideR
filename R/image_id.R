#' Derive a per-pixel image identity from a peak matrix
#'
#' @description
#' Internal helper returning a factor with one entry per pixel (row of
#' `pkm$intensity`) that identifies which image/run each pixel belongs to, so
#' that peak matrices coming from different input paths share one
#' representation.
#'
#' The identity is taken, in order of preference, from:
#' \enumerate{
#'   \item `pkm$image_id`, as added by [cardinal_to_peakmatrix()].
#'   \item `pkm$numPixels` and `pkm$names`, as carried natively by rMSI2
#'     `rMSIprocPeakMatrix` objects with more than one image: pixels are
#'     concatenated in contiguous row-blocks following that order.
#' }
#' If neither describes more than one image, every pixel is assigned to a
#' single image (named after `pkm$names[1]` when available).
#'
#' @param pkm Peak matrix list containing at least `intensity`.
#'
#' @return A factor of length `nrow(pkm$intensity)`, with levels in order of
#'   first appearance.
#'
#' @keywords internal
derive_image_id_from_pkm <- function(pkm) {
  if (!is.list(pkm) || !is.matrix(pkm$intensity)) {
    stop("`pkm` must be a list containing a numeric `intensity` matrix.", call. = FALSE)
  }

  n_pix <- nrow(pkm$intensity)

  if (!is.null(pkm$image_id)) {
    if (length(pkm$image_id) != n_pix) {
      stop(
        "Length of `pkm$image_id` must match the number of rows in `pkm$intensity`.",
        call. = FALSE
      )
    }
    ids <- as.character(pkm$image_id)
    return(factor(ids, levels = unique(ids)))
  }

  n_pixels <- pkm$numPixels
  img_names <- as.character(pkm$names)

  if (length(n_pixels) > 1L) {
    if (sum(n_pixels) != n_pix) {
      stop(
        "Sum of `pkm$numPixels` must match the number of rows in `pkm$intensity`.",
        call. = FALSE
      )
    }
    if (length(img_names) != length(n_pixels)) {
      img_names <- paste0("image_", seq_along(n_pixels))
    }
    img_names <- make.unique(img_names)
    return(factor(rep(img_names, times = n_pixels), levels = img_names))
  }

  single_name <- if (length(img_names) >= 1L && !is.na(img_names[1]) && nzchar(img_names[1])) {
    img_names[1]
  } else {
    "image_1"
  }
  factor(rep(single_name, n_pix))
}
