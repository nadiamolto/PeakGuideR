
# ============================================================================================
# Isotope detection module:Only those peaks that have a high score will pass to the next step
# ============================================================================================
#
# This function computes R^2 colocalization between M+0 and its isotopic
# candidates (M+1..M+K) across pixels. It supports:
#   - 13C series (Δ ≈ 1.003355 / z) for k = 1..max_iso (maximum number the isotopes to look for)
#   - optional "heavy" ~2 Da M+2 candidates (Cl/Br/S) at Δ ≈ 2 / z
# They will be considerated depending on colocalization and mass dif)
# It validates data, filters paired zeros, and combines the per-step R^2
# into a final morphology score using min/mean/weighted schemes.
#
# INPUT:
#   peakmatrix: pkm$mas (list), pkm$intensity (matrix where pixels=rows)
#
# KEY PARAMETERS:
#   max_iso            how many isotope steps to consider (1 = M+1 only; 2 = M+1 & M+2; ...)
#   charges            vector of charge states to try (e.g., c(1,2))
#   include_heavy_M2   if TRUE, for k==2 also test ~2Da shifts (Cl/Br/S) along with 13C2
#   mode               candidate search mode: "ppm" (centroided) or "dp" (profile/binned)
#   tol_ppm / tol_dp   tolerance window in ppm or data points.
#   require_contiguity if TRUE, only accept M+2 if M+1 was valid
#   combine            how to combine per-step R^2 into a single score: "min" | "mean" | "weighted"
#   weights            numeric vector of length max_iso (only used if combine == "weighted")
#
# OUTPUT:
#   A data.frame with (per M0 feature):
#     idx_M0, mz_M0
#     idx_Mk, mz_Mk, r2_Mk  for k = 1..max_iso
#     score_morph           (combined morphology score across k)
#
# Notes:
#   - All indices refer to the (possibly re-ordered) local mass vector (we sort by m/z).
#   - You can map back to original order if needed (see comment near the end).
#

# ---- Internal helpers --------------------------------------

.check_peakmatrix <- function(pm) {
  stopifnot(is.list(pm), !is.null(pm$mass), !is.null(pm$intensity))
  stopifnot(is.numeric(pm$mass), is.matrix(pm$intensity))
  if (length(pm$mass) != ncol(pm$intensity)) {
    stop("length(pm$mass) must equal ncol(pm$intensity).", call. = FALSE)
  }
}

.ppm_tol <- function(target_mz, ppm) target_mz * ppm * 1e-6 # it calculates ppm in Da

# Find closest candidate by ppm window
######## IMPORTANT: It is recommended to use it for centered data or without reliable bin size 
.find_candidate_ppm <- function(masses, target, ppm) {
  tol <- .ppm_tol(target, ppm) #calculate tolerance in Da depending on the ppm introduced
  diffs <- abs(masses - target) #absolutedistance 
  ok <- which(diffs <= tol) #columns where the value is into the mass+-tol
  if (!length(ok)) return(NA_integer_)
  ok[which.min(diffs[ok])] #Among all the valid data, the nearest is saved
}

# Find closest candidate by "data point" (dp)
########## IMPORTANT: Use it when the matrix is regularly binned (spaces between columns are constant)
.find_candidate_dp <- function(masses, target, tol_dp) {
  # Nearest index in sorted 'masses' (they are sorted previously)
  k <- findInterval(target, masses)
  # Candidate set: nearest ± tol_dp steps (clip to [1, length])
  cand <- seq.int(max(1L, k - tol_dp), min(length(masses), k + tol_dp)) #window
  if (!length(cand)) return(NA_integer_)
  diffs <- abs(masses[cand] - target)
  cand[which.min(diffs)] #nearest candidate 
}


# Considering bin size. We do not assume a global regularity in feature width
#The best result will be the nearest feature within +- N bins * binsize (Da)
#### IMPORTANT: THis is the best method for binned matrix with no regularity
.find_candidate_binsize <- function(masses, binsizes, target, N_bins = 3L) {
  idx_nearest <- which.min(abs(masses - target))
  tol <- N_bins * binsizes[idx_nearest] #Tolerance depends on the bin size 
  diffs <- abs(masses - target)
  ok <- which(diffs <= tol)
  if (!length(ok)) return(NA_integer_)
  ok[which.min(diffs[ok])]
}


# Compute R^2 between y (candidate) and x (base M+0) with robust validations:
# paired-zero filtering (drop pixels where x==0 & y==0)
# require minimum pixel count and non-zero variance
# choose intercept vs. through-origin regression
.r2_colocalization <- function(x, y, min_pixels = 20L, use_intercept = TRUE) { #at least 20 pixels to perform it
  # Paired-zero filtering: drop pixels where BOTH channels are zero
  keep <- !(x == 0 & y == 0)
  # Also enforce finite values
  keep <- keep & is.finite(x) & is.finite(y)
  if (sum(keep) < max(3L, min_pixels)) return(NA_real_) #If the number of pixel is lower than 3 or
  xk <- x[keep]; yk <- y[keep]
  if (var(xk) == 0 || var(yk) == 0) return(NA_real_)
  
  df <- data.frame(y = yk, x = xk)
  fit <- if (use_intercept) stats::lm(y ~ x, data = df) else stats::lm(y ~ 0 + x, data = df)
  summary(fit)$r.squared
}

# Build list of isotope mass deltas (in Da) to test for a given k and charge(s).
# For 13C series: k * 1.003355 / z
# For "heavy" M+2 (k==2 only, if include_heavy_M2): add ~2 Da deltas / z for Cl/Br/S
.isotope_deltas_for_k <- function(k, charges, include_heavy_M2) {
  # 13C shift per step (Daltons)
  d13C <- 1.003355
  base <- as.numeric(outer(k * d13C, 1 / charges))  # k * 1.003355 / z
  
  if (k == 2L && isTRUE(include_heavy_M2)) {
    # Add "heavy" ~2Da alternatives (divide by z as well)
    heavy_raw <- c(1.99705, 1.99795, 1.99580)  # 37Cl, 81Br, 34S
    heavy <- as.numeric(outer(heavy_raw, 1 / charges))  # Δ / z
    deltas <- c(base, heavy)
  } else {
    deltas <- base
  }
  sort(unique(deltas))
}

# ---- Main function -----------------------------------------------------------

#' Isotope morphology score using binSize-aware tolerance (preferred for binned data)
#'
#' @param pm list with fields:
#'   - mass: numeric vector of length p (feature m/z centers).
#'   - intensity: numeric matrix [Pixels x p].
#'   - binSize: (optional) numeric vector of length p (bin width in Da per feature).
#' @param prefer_mode character: "binsize" (default), "dp", or "ppm".
#'   If "binsize" is chosen but binSize is missing, it falls back to "dp" (then "ppm").
#' @param N_bins integer: multiplier for binSize tolerance (± N_bins * binSize) when using binsize mode.
#' @param tol_dp integer: ± data points window (when using dp mode).
#' @param tol_ppm numeric: ppm window (when using ppm mode).
#' @param max_iso integer: how many isotope steps (k) to test (1 = M+1; 2 = M+1 & M+2; ...).
#' @param charges integer vector: charge states to consider (e.g., c(1L,2L)).
#' @param include_heavy_M2 logical: if TRUE, also test ~2 Da M+2 shifts (Cl/Br/S).
#' @param require_contiguity logical: if TRUE, only accept M+2 if M+1 was valid (non-NA R^2).
#' @param min_pixels integer: minimum pixels used in regression after paired-zero filtering.
#' @param use_intercept logical: TRUE → lm(y ~ x); FALSE → lm(y ~ 0 + x).
#' @param combine "min"|"mean"|"weighted": how to combine per-step R^2 into final morphology score.
#' @param weights numeric: weights for k=1..max_iso if combine="weighted".
#'
#' @return data.frame with per-M0 results:
#'   idx_M0, mz_M0, score_morph, and for k=1..max_iso: idx_Mk, mz_Mk, r2_Mk
#' @export

iso_morphology_series_binsafe <- function(
    pm,
    prefer_mode = c("binsize", "dp", "ppm"),
    N_bins      = 3L,
    tol_dp      = 3L,
    tol_ppm     = 5,
    max_iso     = 2L,
    charges     = 1L,
    include_heavy_M2   = TRUE,
    require_contiguity = TRUE,
    min_pixels  = 20L,
    use_intercept = TRUE,
    combine = c("min", "mean", "weighted"),
    weights = NULL
) {
  .check_peakmatrix(pm)
  prefer_mode <- match.arg(prefer_mode)
  combine <- match.arg(combine)
  
  mass <- as.numeric(pm$mass)
  Imat <- pm$intensity
  p <-ncol(Imat)
  
  # Sort by m/z for efficient neighborhood search
  ord <- order(mass)
  mass <- mass[ord]
  Imat <- Imat[, ord, drop = FALSE]
  binSize <- if (!is.null(pm$binSize)) as.numeric(pm$binSize)[ord] else NULL
 
  
  # Decide finder based on availability and preference
  if (prefer_mode == "binsize" && !is.null(binSize)) {
    finder <- function(target) .find_candidate_binsize(mass, binSize, target, N_bins = N_bins)
    err_metric <- function(target, found) abs(mass[found] - target) # Da
  } else if (prefer_mode %in% c("binsize", "dp")) {
    # fallback to dp if binSize is not available
    finder <- function(target) .find_candidate_dp(mass, target, tol_dp)
    err_metric <- function(target, found) abs(mass[found] - target) # Da
  } else {
    # ppm mode
    finder <- function(target) .find_candidate_ppm(mass, target, tol_ppm)
    err_metric <- function(target, found) abs(1e6 * (mass[found] - target) / target) # ppm
  }
  
  out <- vector("list", p)
  
  for (i in seq_len(p)) {
    mz0 <- mass[i]
    I0  <- Imat[, i]
    
    idx_k <- rep(NA_integer_, max_iso)
    mz_k  <- rep(NA_real_,     max_iso)
    r2_k  <- rep(NA_real_,     max_iso)
    
    prev_valid <- TRUE
    
    for (k in seq_len(max_iso)) {
      if (require_contiguity && !prev_valid) break
      
      deltas <- .isotope_deltas_for_k(k, charges, include_heavy_M2)
      best <- NULL
      
      for (d in deltas) {
        target <- mz0 + d
        j <- finder(target)
        if (is.na(j)) next
        
        r2 <- .r2_colocalization(I0, Imat[, j], min_pixels = min_pixels, use_intercept = use_intercept)
        err <- err_metric(target, j)
        
        cand <- list(j = j, mz = mass[j], r2 = r2, err = err)
        if (is.null(best)) {
          best <- cand
        } else {
          # Prefer a valid r2; if tie (both NA or both valid), pick the smaller error
          better <- if (is.na(best$r2) && !is.na(cand$r2)) TRUE else
            if (!is.na(best$r2) && is.na(cand$r2)) FALSE else
              (cand$err < best$err)
          if (better) best <- cand
        }
      }
      
      if (!is.null(best)) {
        idx_k[k] <- best$j
        mz_k[k]  <- best$mz
        r2_k[k]  <- best$r2
        prev_valid <- !is.na(best$r2)
      } else {
        prev_valid <- FALSE
      }
    } # k
    
    # Combine per-step R^2
    r2_valid <- r2_k[!is.na(r2_k)]
    score_morph <- if (!length(r2_valid)) NA_real_ else switch(
      combine,
      "min" = min(r2_valid),
      "mean" = mean(r2_valid),
      "weighted" = {
        if (is.null(weights)) {
          # default: emphasize lower k (M+1 > M+2 > M+3)
          w <- rev(seq_len(max_iso))
        } else {
          stopifnot(length(weights) >= max_iso)
          w <- weights[seq_len(max_iso)]
        }
        w_use <- w[!is.na(r2_k)]
        w_use <- w_use / sum(w_use)
        sum(r2_valid * w_use)
      }
    )
    
    out[[i]] <- data.frame(
      idx_M0 = i, mz_M0 = mz0,
      score_morph = score_morph,
      idx_M1 = if (max_iso >= 1) idx_k[1] else NA_integer_,
      mz_M1  = if (max_iso >= 1) mz_k[1]  else NA_real_,
      r2_M1  = if (max_iso >= 1) r2_k[1]  else NA_real_,
      idx_M2 = if (max_iso >= 2) idx_k[2] else NA_integer_,
      mz_M2  = if (max_iso >= 2) mz_k[2]  else NA_real_,
      r2_M2  = if (max_iso >= 2) r2_k[2]  else NA_real_,
      idx_M3 = if (max_iso >= 3) idx_k[3] else NA_integer_,
      mz_M3  = if (max_iso >= 3) mz_k[3]  else NA_real_,
      r2_M3  = if (max_iso >= 3) r2_k[3]  else NA_real_
    )
  }
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}