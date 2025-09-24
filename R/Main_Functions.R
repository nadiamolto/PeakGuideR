
# Internal helper: comvert MSImagingExperiment object to a peakmatrix object (keeping the treacability to include the
# metadata in the final output)

cardinal_to_rMSIprocPeakMatrix <- function(
    x,
    value = NULL,        # <-- ahora es NULL por defecto
    dataset_name = NULL,
    snr = NULL,
    area = NULL
) {
  if (!inherits(x, "MSImagingExperiment")) {
    stop("`x` must be a MSImagingExperiment (Cardinal ≥ 3.x).", call. = FALSE)
  }
  
  # --- pick layer automatically if value is NULL ---
  img <- Cardinal::imageData(x)
  avail <- names(img)
  if (is.null(value)) {
    if ("intensity" %in% avail) {
      value <- "intensity"
    } else if (length(avail) >= 1L) {
      value <- avail[1]
      message(sprintf("`value` not provided. Using first available layer: '%s'", value))
    } else {
      stop("imageData(x) has no layers.", call. = FALSE)
    }
  }
  if (!value %in% avail) {
    stop(sprintf("Layer '%s' not found in imageData(x). Available: %s",
                 value, paste(avail, collapse = ", ")), call. = FALSE)
  }
  
  # --- intensity layer: features x spectra -> transpose to pixels x features ---
  intens_layer <- img[[value]]
  if (is.matrix(intens_layer)) {
    M <- intens_layer
  } else if (methods::is(intens_layer, "DelayedArray") || methods::is(intens_layer, "ImageArray")) {
    M <- as.matrix(intens_layer)
  } else {
    M <- try(as.matrix(intens_layer), silent = TRUE)
    if (inherits(M, "try-error")) stop("Could not coerce selected layer to matrix.", call. = FALSE)
  }
  intensity <- t(M)  # pixels x features
  
  # --- mass (m/z) vector ---
  mass <- tryCatch({
    fd <- as.data.frame(Cardinal::featureData(x))
    if ("mz" %in% names(fd)) as.numeric(fd$mz) else as.numeric(Cardinal::mz(x))
  }, error = function(e) {
    as.numeric(Cardinal::mz(x))
  })
  if (length(mass) != ncol(intensity)) {
    warning("Length of 'mass' does not match ncol(intensity). Creating a simple sequence.")
    mass <- seq_len(ncol(intensity))
  }
  
  # --- binSize heuristic (local spacing / 2; endpoints use neighbor spacing) ---
  binSize <- rep(NA_real_, length(mass))
  if (length(mass) >= 2) {
    dm <- diff(mass)
    if (length(mass) == 2) {
      binSize[] <- dm[1] / 2
    } else {
      left  <- c(dm[1], dm)
      right <- c(dm, dm[length(dm)])
      binSize <- pmin(left, right) / 2
    }
  }
  
  # --- pixel positions (x, y) ---
  pd <- tryCatch(as.data.frame(Cardinal::pixelData(x)),
                 error = function(e) stop("Could not access pixelData(x).", call. = FALSE))
  if (!("x" %in% names(pd) && "y" %in% names(pd))) {
    stop("pixelData(x) must contain 'x' and 'y' columns.", call. = FALSE)
  }
  pos <- cbind(x = as.numeric(pd$x), y = as.numeric(pd$y))
  rownames(pos) <- NULL
  posMotors <- pos
  
  # --- per-pixel normalizations ---
  TIC <- rowSums(intensity, na.rm = TRUE)
  RMS <- sqrt(rowMeans(intensity^2, na.rm = TRUE))
  MAX <- apply(intensity, 1L, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  normalizations <- data.frame(TIC = TIC, RMS = RMS, MAX = MAX, row.names = NULL)
  
  # --- SNR & area (zeros by default unless provided) ---
  nPix <- nrow(intensity); nFeat <- ncol(intensity)
  if (is.null(snr))  snr  <- matrix(0, nrow = nPix, ncol = nFeat)
  if (is.null(area)) area <- matrix(0, nrow = nPix, ncol = nFeat)
  
  # --- dataset name & uuid ---
  if (is.null(dataset_name)) {
    rn <- tryCatch(Cardinal::runNames(x), error = function(e) NULL)
    dataset_name <- if (!is.null(rn) && length(rn) >= 1 && nzchar(rn[1])) rn[1] else "cardinal_import"
  }
  uuid <- .make_uuid_hex()
  
  # --- assemble rMSIprocPeakMatrix ---
  pkm_rms <- list(
    mass           = mass,
    binSize        = binSize,
    intensity      = intensity,
    SNR            = snr,
    area           = area,
    normalizations = normalizations,
    pos            = pos,
    numPixels      = as.integer(nPix),
    names          = as.character(dataset_name),
    uuid           = uuid,
    posMotors      = posMotors
  )
  colnames(pkm_rms$pos)       <- c("x", "y")
  colnames(pkm_rms$posMotors) <- c("x", "y")
  class(pkm_rms) <- "rMSIprocPeakMatrix"
  pkm_rms
}