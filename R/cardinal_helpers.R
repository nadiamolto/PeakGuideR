#' Check whether an object is a supported Cardinal object
#'
#' @description
#' Internal helper used to detect Cardinal MSI objects before conversion to a
#' PeakGuideR peak matrix.
#'
#' @param x An R object.
#'
#' @return Logical value.
#'
#' @keywords internal
is_cardinal_object <- function(x) {
  supported_classes <- c(
    "MSImagingExperiment",
    "MSContinuousImagingExperiment",
    "MSProcessedImagingExperiment"
  )

  any(vapply(
    supported_classes,
    function(cls) {
      inherits(x, cls) || methods::is(x, cls)
    },
    logical(1)
  ))
}


#' Convert a Cardinal MSImagingExperiment object to a peak matrix
#'
#' @description
#' Converts a supported Cardinal MSI object into a PeakGuideR peak matrix with
#' an rMSIprocPeakMatrix-like structure.
#'
#' @param x A Cardinal MSImagingExperiment object.
#' @param value Name of the imageData layer to extract. If NULL, "intensity" is used.
#' @param dataset_name Optional dataset name.
#' @param snr Optional SNR matrix.
#' @param area Optional area matrix.
#'
#' @return An object with rMSIprocPeakMatrix-like structure.
#'
#' @examples
#' \dontrun{
#' # `cardinal_msi_data` is a Cardinal MSImagingExperiment object.
#' pkm <- cardinal_to_peakmatrix(cardinal_msi_data)
#'
#' res <- run_peakguider_workflow(
#'   pkm = pkm,
#'   ion_mode = "pos",
#'   matrix = "HCCA"
#' )
#' }
#'
#' @export
cardinal_to_peakmatrix <- function(
    x,
    value = NULL,
    dataset_name = NULL,
    snr = NULL,
    area = NULL
) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop(
      "Package 'Cardinal' is required to use `cardinal_to_peakmatrix()`. ",
      "Please install it with: BiocManager::install('Cardinal')",
      call. = FALSE
    )
  }

  if (!is_cardinal_object(x)) {
    stop(
      "`x` must be a supported Cardinal MSImagingExperiment object.",
      call. = FALSE
    )
  }

  img <- Cardinal::imageData(x)
  avail <- names(img)

  if (is.null(value)) {
    if (!("intensity" %in% avail)) {
      stop("imageData(x) does not contain an 'intensity' layer.", call. = FALSE)
    }
    value <- "intensity"
  }

  if (!value %in% avail) {
    stop(sprintf(
      "Layer '%s' not found in imageData(x). Available layers: %s",
      value,
      paste(avail, collapse = ", ")
    ), call. = FALSE)
  }

  intens_layer <- img[[value]]

  if (is.matrix(intens_layer)) {
    M <- intens_layer
  } else if (
    methods::is(intens_layer, "DelayedArray") ||
    methods::is(intens_layer, "ImageArray")
  ) {
    M <- as.matrix(intens_layer)
  } else {
    M <- try(as.matrix(intens_layer), silent = TRUE)
    if (inherits(M, "try-error")) {
      # Some S4 classes used by Cardinal's "processed"/peak-picked
      # imaging objects (e.g. matter::sparse_mat) register an as.matrix()
      # S4 method, but standardGeneric() dispatch for that method is only
      # visible to callers whose namespace imports the defining package.
      # Cardinal/matter are optional dependencies (Suggests, guarded by
      # requireNamespace() above), so PeakGuideR's own namespace never
      # imports them and the implicit dispatch above can fail even though
      # the method exists and is registered globally. Fall back to an
      # explicit lookup in the S4 method table, which does not depend on
      # the caller's namespace imports.
      m <- tryCatch(
        methods::selectMethod("as.matrix", class(intens_layer), optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(m)) {
        M <- try(m(intens_layer), silent = TRUE)
      }
    }
    if (inherits(M, "try-error")) {
      stop("Could not coerce selected imageData layer to matrix.", call. = FALSE)
    }
  }

  intensity <- t(M)

  mass <- tryCatch({
    fd <- as.data.frame(Cardinal::featureData(x))
    if ("mz" %in% names(fd)) {
      as.numeric(fd$mz)
    } else {
      as.numeric(Cardinal::mz(x))
    }
  }, error = function(e) {
    as.numeric(Cardinal::mz(x))
  })

  if (length(mass) != ncol(intensity)) {
    warning(
      "Length of 'mass' does not match ncol(intensity). ",
      "Creating a simple sequence instead."
    )
    mass <- seq_len(ncol(intensity))
  }

  binSize <- rep(NA_real_, length(mass))

  if (length(mass) >= 2) {
    dm <- diff(mass)

    if (length(mass) == 2) {
      binSize[] <- dm[1] / 2
    } else {
      left <- c(dm[1], dm)
      right <- c(dm, dm[length(dm)])
      binSize <- pmin(left, right) / 2
    }
  }

  pd <- tryCatch(
    as.data.frame(Cardinal::pixelData(x)),
    error = function(e) {
      stop("Could not access pixelData(x).", call. = FALSE)
    }
  )

  if (!all(c("x", "y") %in% names(pd))) {
    stop("pixelData(x) must contain 'x' and 'y' columns.", call. = FALSE)
  }

  pos <- cbind(
    x = as.numeric(pd$x),
    y = as.numeric(pd$y)
  )

  rownames(pos) <- NULL
  posMotors <- pos

  TIC <- rowSums(intensity, na.rm = TRUE)
  RMS <- sqrt(rowMeans(intensity^2, na.rm = TRUE))

  MAX <- apply(
    intensity,
    1L,
    function(v) {
      if (all(is.na(v))) {
        NA_real_
      } else {
        max(v, na.rm = TRUE)
      }
    }
  )

  normalizations <- data.frame(
    TIC = TIC,
    RMS = RMS,
    MAX = MAX,
    row.names = NULL
  )

  nPix <- nrow(intensity)
  nFeat <- ncol(intensity)

  if (is.null(snr)) {
    snr <- matrix(0, nrow = nPix, ncol = nFeat)
  }

  if (is.null(area)) {
    area <- matrix(0, nrow = nPix, ncol = nFeat)
  }

  if (is.null(dataset_name)) {
    rn <- tryCatch(Cardinal::runNames(x), error = function(e) NULL)

    dataset_name <- if (!is.null(rn) && length(rn) >= 1 && nzchar(rn[1])) {
      rn[1]
    } else {
      "cardinal_import"
    }
  }

  pkm_rms <- list(
    mass = mass,
    binSize = binSize,
    intensity = intensity,
    SNR = snr,
    area = area,
    normalizations = normalizations,
    pos = pos,
    numPixels = as.integer(nPix),
    names = as.character(dataset_name),
    posMotors = posMotors
  )

  colnames(pkm_rms$pos) <- c("x", "y")
  colnames(pkm_rms$posMotors) <- c("x", "y")

  class(pkm_rms) <- "rMSIprocPeakMatrix"

  pkm_rms
}
