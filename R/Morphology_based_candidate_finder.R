#' Isotopic morphology candidates: based on co-localization (z = 1 fixed)
#'
#' @description
#' Detects potential isotopic relationships between peaks in a Mass Spectrometry
#' Imaging (MSI) peak matrix based on expected mass differences and spatial
#' colocalization (correlation). It computes correlation-like morphology metrics (Pearson,
#' cosine, or Spearman), optionally evaluates spatial consistency across tiles,
#' and returns all isotopic pairs that satisfy mass tolerance and a minimum
#' similarity threshold.
#'
#' @details
#' Charge state is **fixed to z = 1** throughout the method (no user control).
#'
#' For each monoisotopic feature (M+0), the function searches isotopic
#' candidates at theoretical mass differences for z = 1:
#'
#' - **13C series:** M+1 = 1.003355, M+2 = 2 × 1.003355  
#' - **Other isotopes (M+X level):**  
#'   34S = 1.99580, 37Cl = 1.99705, 81Br = 1.99795, 18O = 2.004245
#'
#' **Tolerance modes:**
#' - `"ppm"` (default): find candidates within ±`tol_ppm` ppm around the theoretical mass.
#' - `"dp"`: find candidates within ±`tol_dp` data points in the ordered mass axis.
#'
#' **Preprocessing steps per intensity pair:**
#' 1. Remove paired zeros (`x==0 & y==0`).
#' 2. Exclude pixels where both intensities lie below the pooled `min_quantile`.
#' 3. Optionally clip negatives (`clip_negatives=TRUE`). 
#' 4. Apply intensity transform (`none`, `log1p`, or `zscore`). Scoring metrics should be considered before setting 
#'    the transformation method: zscore is recommended for pearson correlation, log1p or none for cosine similarity
#'    and none for spearman.
#'    Normalization can be applied (e.g. TIC) before it but should be considered regarding result interpretation.
#'
#' **Scoring metrics (`method`):**
#' - `"pearson"`: squared Pearson correlation (R²), [0,1].
#' - `"cosine"`: cosine similarity.
#' - `"spearman"`: squared Spearman rank correlation (negatives set to 0).
#'
#' **Tile mode:** if `use_tiles=TRUE` and pixel coordinates are present (`pm$pos$x,y`),
#' the sample is split into approximately `sqrt(tiles)×sqrt(tiles)` subregions.
#' For each tile, the same correlation is computed, summarizing via:
#' - `"median"`: median of tile scores,
#' - `"p25"`: 25th percentile (robust lower-bound),
#' - `"pass_rate"`: fraction of tiles with score ≥ `tile_threshold`.
#'
#' A consistency factor is computed as:
#' \deqn{ tile\_consistency = \max(0, 1 - sd(tile\_scores)) }
#'
#' and blended with the global score:
#' \deqn{ score\_final = \alpha \cdot score\_global + (1-\alpha) \cdot tile\_consistency }
#' Alpha will be set by the user giving weights to global and tile analysis.
#'
#' **Output policy:**
#' Returns *all* isotopic candidate pairs (M+0... M+k) that fall within the mass
#' tolerance and have `score_final ≥ min_score_keep`.
#' A dataframe with as many rows as isotopes are found.
#'
#' @param pm A list with:
#'   - `mass`: numeric vector of m/z values.
#'   - `intensity`: matrix [pixels × features].
#'   - optional `pos`: data.frame/matrix with `x` and `y` coordinates for pixels
#'     (without it, tile strategy is disabled with an info message).
#'     Position information is recommended to get more reliable results.
#' @param prefer_mode `"ppm"` (default) or `"dp"`. Tolerance mode.
#' @param tol_ppm Numeric, mass window (±ppm). Default `5`.
#' @param tol_dp Integer, window in indices if `prefer_mode="dp"`. Default `3L`.
#' @param method Character, morphology metric: `"pearson"` (default), `"cosine"`, or `"spearman"`.
#' @param use_intercept Logical, reserved for compatibility (currently unused). Default `FALSE`.
#' @param transform One of `"none"`, `"log1p"`, `"zscore"`. Default `"log1p"`.
#' @param min_quantile Numeric in [0,1], pooled lower quantile filter. Default `0.01`.
#' @param clip_negatives Logical, clip negative values before transform. Default `TRUE`.
#' @param use_tiles Logical, enable spatial tile evaluation if `pm$pos` exists. Default `TRUE`.
#' @param tiles Integer, approximate number of subregions. Default `9L`.
#' @param tile_min_pixels Integer, minimum pixels per tile. Default `50L`.
#' @param tile_blend One of `"median"`, `"p25"`, or `"pass_rate"`. Default `"median"`.
#' @param tile_alpha Numeric in [0,1], blending between global and tile score. Default `0.8`.
#' @param tile_threshold Numeric in [0,1], threshold for `"pass_rate"`. Default `0.8`.
#' @param min_score_keep Numeric in [0,1], score threshold. Default `0.2`.
#'
#' @return A `data.frame` with:
#' \itemize{
#'   \item `idx_M0`, `mz_M0`: anchor index and m/z.
#'   \item `iso_type`: isotope type (`C13_M1`, `C13_2`, `S34`, `Cl37`, `Br81`, `O18`).
#'   \item `k`: isotope level (1 or 2).
#'   \item `z`: charge state (always 1).
#'   \item `idx_cand`, `mz_cand`: candidate index and m/z.
#'   \item `score_global`, `tile_summary`, `tile_sd`, `tile_consistency`, `score_final`.
#'   \item `mass_err_da`, `mass_err_ppm`: candidate error vs theoretical mass.
#' }
#'
#' @examples
#' \dontrun{
#' res <- iso_morphology_candidates(
#'   pm,
#'   prefer_mode   = "ppm",
#'   tol_ppm       = 5,
#'   method        = "pearson",
#'   transform     = "log1p",
#'   use_tiles     = TRUE,
#'   tiles         = 9L,
#'   tile_blend    = "median",
#'   tile_alpha    = 0.8,
#'   min_score_keep= 0.2
#' )
#' head(res)
#' }
#' @export
iso_morphology_candidates <- function(
    pm,
    prefer_mode    = c("ppm","dp"),
    tol_ppm        = 5,
    tol_dp         = 3L,
    method         = c("pearson","cosine","spearman"),
    use_intercept  = FALSE,
    transform      = c("none","log1p","zscore"),
    min_quantile   = 0.01,
    clip_negatives = TRUE,
    use_tiles      = TRUE,
    tiles          = 9L,
    tile_min_pixels= 50L,
    tile_blend     = c("median","p25","pass_rate"),
    tile_alpha     = 0.8,
    tile_threshold = 0.8,
    min_score_keep = 0.2
) {
  # ---- Argument normalization ----
  prefer_mode <- match.arg(prefer_mode)
  method      <- match.arg(method)
  transform   <- match.arg(transform)
  tile_blend  <- match.arg(tile_blend)
  
  # ---- Validate input ----
  stopifnot(is.list(pm), !is.null(pm$mass), !is.null(pm$intensity))
  mass <- as.numeric(pm$mass)
  Imat <- pm$intensity
  stopifnot(is.numeric(mass), is.matrix(Imat), length(mass) == ncol(Imat))
  
  # ---- Handle coordinates ----
  pos_xy <- if (!is.null(pm$pos) && all(c("x","y") %in% colnames(pm$pos))) pm$pos else NULL
  if (isTRUE(use_tiles) && is.null(pos_xy)) {
    message("[INFO] No 'x','y' coordinates found in pm$pos. Tile analysis disabled.")
    use_tiles <- FALSE
  }
  
  # ---- Sort by m/z ----
  ord   <- order(mass) #From lower to higher
  mass  <- mass[ord]  #Mass reordering
  Imat  <- Imat[, ord, drop=FALSE] #Matrix reordering including intensities
  p     <- ncol(Imat) #Number of masses
  
  # ---- Isotopic deltas for z = 1 (fixed) ----
  iso_table <- (function() {
    d13C  <- 1.003355
    heavy <- c(S34=1.99580, Cl37=1.99705, Br81=1.99795, O18=2.004245)
    data.frame(
      iso_type = c("C13_M1", "C13_2", names(heavy)),
      k        = c(1L, 2L, rep(2L, length(heavy))),
      z        = 1L,
      delta    = c(d13C, 2*d13C, unname(heavy)),
      stringsAsFactors = FALSE
    )
  })()
  
  # ---- Candidate finder ----
  get_candidates <- if (prefer_mode == "ppm") {
    function(target) {
      tol <- target * tol_ppm * 1e-6
      which(abs(mass - target) <= tol) #Lower than the tolerance 
    }
  } else {
    function(target) {
      k <- findInterval(target, mass) #In which interval the target is 
      seq.int(max(1L, k - tol_dp), min(p, k + tol_dp)) #index sequence generator into the tolerance window
    }
  }
  
  # ---- Preprocessing ----
  preprocess_xy <- function(x, y) {
    keep <- is.finite(x) & is.finite(y) & !(x == 0 & y == 0) #boolean
    if (min_quantile > 0 && any(keep)) { #quantil filter
      pool <- c(x[keep], y[keep]) #keep if true
      cutoff <- stats::quantile(pool, probs=min_quantile, na.rm=TRUE, type=7) #lowest quantil using package stats
      keep <- keep & !(x <= cutoff & y <= cutoff) #remove pixels under cutoff
    }
    if (sum(keep) < 3L) return(NULL) #more than 3
    xk <- x[keep]; yk <- y[keep]
    if (clip_negatives) { xk <- pmax(xk, 0); yk <- pmax(yk, 0) }
    tf <- switch(transform,
                 "none" = identity,
                 "log1p" = log1p,
                 "zscore" = function(v) {
                   sdv <- stats::sd(v)
                   if (is.na(sdv) || sdv == 0) rep(0, length(v)) else (v - mean(v))/sdv
                 })
    list(x=tf(xk), y=tf(yk))
  }
  
  # ---- Scoring ----
  score_core <- function(x, y) {
    if (length(x) != length(y) || length(x) < 3L) return(NA_real_)
    if (stats::var(x)==0 || stats::var(y)==0) return(NA_real_)
    switch(method,
           "pearson" = {
             r <- suppressWarnings(stats::cor(x, y, method="pearson"))
             if (is.na(r)) NA_real_ else max(0, min(1, r*r))
           },
           "cosine" = {
             num <- sum(x*y); den <- sqrt(sum(x^2))*sqrt(sum(y^2))
             if (den==0) NA_real_ else max(0, min(1, num/den))
           },
           "spearman" = {
             r <- suppressWarnings(stats::cor(x, y, method="spearman"))
             if (is.na(r) || r <= 0) 0 else max(0, min(1, r*r))
           })
  }
  
  # ---- Tile scoring ----
  split_into_tiles <- function(pos_xy, n_tiles) {
    s <- max(1L, round(sqrt(n_tiles)))
    xbreaks <- unique(stats::quantile(pos_xy[, "x"], seq(0,1,length.out=s+1), na.rm=TRUE, type=7))
    ybreaks <- unique(stats::quantile(pos_xy[, "y"], seq(0,1,length.out=s+1), na.rm=TRUE, type=7))
    if (length(xbreaks) < 2L || length(ybreaks) < 2L) {
      return(rep("1_1", nrow(pos_xy)))
    }
    tx <- cut(pos_xy[, "x"], breaks = xbreaks, include.lowest = TRUE, labels = FALSE)
    ty <- cut(pos_xy[, "y"], breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    paste(tx, ty, sep="_")
  }
  
  tile_scores <- function(I0, Ik) {
    ids <- split_into_tiles(pos_xy, tiles)
    uids <- unique(ids)
    scores <- numeric(length(uids))
    for (kk in seq_along(uids)) {
      idx <- which(ids == uids[kk])
      if (length(idx) < tile_min_pixels) { scores[kk] <- NA_real_; next }
      pp <- preprocess_xy(I0[idx], Ik[idx])
      if (is.null(pp)) { scores[kk] <- NA_real_; next }
      scores[kk] <- score_core(pp$x, pp$y)
    }
    list(
      median     = stats::median(scores, na.rm=TRUE),
      p25        = stats::quantile(scores, 0.25, na.rm=TRUE, type=7),
      pass_rate  = mean(scores >= tile_threshold, na.rm=TRUE),
      sd         = stats::sd(scores, na.rm=TRUE)
    )
  }
  
  # ---- Iterate over anchors ----
  rows <- list(); rr <- 0L
  for (i in seq_len(p)) {
    mz0 <- mass[i]; I0 <- Imat[, i]
    for (h in seq_len(nrow(iso_table))) {
      iso_type <- iso_table$iso_type[h]
      kstep    <- iso_table$k[h]
      zfix     <- 1L
      d        <- iso_table$delta[h]
      target   <- mz0 + d
      
      cand_idx <- get_candidates(target)
      cand_idx <- setdiff(cand_idx, i)
      if (!length(cand_idx)) next
      
      for (j in cand_idx) {
        pp <- preprocess_xy(I0, Imat[, j])
        if (is.null(pp)) next
        
        s_global <- score_core(pp$x, pp$y)
        
        tile_summary <- NA_real_; tile_sd <- NA_real_; tile_consistency <- NA_real_
        if (use_tiles && !is.null(pos_xy)) {
          ts <- tile_scores(I0, Imat[, j])
          tile_summary <- switch(tile_blend,
                                 "median"    = ts$median,
                                 "p25"       = ts$p25,
                                 "pass_rate" = ts$pass_rate)
          tile_sd <- ts$sd
          tile_consistency <- if (is.finite(tile_sd)) max(0, 1 - tile_sd) else NA_real_
        }
        
        s_final <- if (is.finite(tile_consistency))
          tile_alpha * s_global + (1 - tile_alpha) * tile_consistency else s_global
        
        if (is.na(s_final) || s_final < min_score_keep) next
        
        err_da  <- mass[j] - target
        err_ppm <- 1e6 * (mass[j] - target) / target
        
        rr <- rr + 1L
        rows[[rr]] <- data.frame(
          idx_M0          = i,
          mz_M0           = mz0,
          iso_type        = iso_type,
          k               = kstep,
          z               = zfix,               # always 1
          idx_cand        = j,
          mz_cand         = mass[j],
          score_global    = s_global,
          tile_summary    = tile_summary,
          tile_sd         = tile_sd,
          tile_consistency= tile_consistency,
          score_final     = s_final,
          mass_err_da     = err_da,
          mass_err_ppm    = err_ppm,
          stringsAsFactors= FALSE
        )
      }
    }
  }
  
  if (!length(rows)) {
    return(data.frame(
      idx_M0=integer(), mz_M0=numeric(), iso_type=character(),
      k=integer(), z=integer(), idx_cand=integer(), mz_cand=numeric(),
      score_global=numeric(), tile_summary=numeric(),
      tile_sd=numeric(), tile_consistency=numeric(),
      score_final=numeric(), mass_err_da=numeric(),
      mass_err_ppm=numeric(), stringsAsFactors=FALSE
    ))
  }
  
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}




