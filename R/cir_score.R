#' CIR score for C13 M+1 validation
#'
#' Filters C13_M1 candidates and validates against theoretical ratios.
#'
#' @param result Isotope detection results (`iso_morphology_candidates()` output)
#' @param pkm Peak matrix list (`mass`, `intensity`)
#' @param min_score_final Minimum `score_final` threshold [0,1]. Default: 0.6
#' @param cir_rel_tol CIR relative tolerance [0,1]. Default: 0.2
#' @param ratio_method `"mean"` (default), `"median"`, or `"sum"` for I_M1/I_M0
#'
#' @return Enhanced data frame:
#' \itemize{
#'   \item `R_obs`: observed I_M1/I_M0 ratio
#'   \item `R_theo`: theoretical CIR from GAM
#'   \item `cir_rel_err`: relative error
#'   \item `is_valid_c13`: cir_rel_err < cir_rel_tol
#'   \item `cir_score`: 1 - cir_rel_err [0,1]
#' }
#' @export
cir_score <- function(result, pkm,
                      min_score_final = 0.6, # 1-extrapolation error
                      cir_rel_tol = 0.2, # interpolation relative error accepted to consider a candidate
                      ratio_method = c("mean", "median", "sum")) {

  ratio_method <- match.arg(ratio_method)

  # Filter ONLY C13_M1 with good morphology
  c13_pairs <- result %>%
    filter(iso_type == "C13_M1",
           score_final >= min_score_final)

  if (nrow(c13_pairs) == 0) {
    message("No C13_M1 candidates with score_final >= ", min_score_final)
    return(c13_pairs)
  }

  # Indices
  idx_M0 <- c13_pairs$idx_M0
  idx_M1 <- c13_pairs$idx_cand

  # Intensities
  I_M0 <- pkm$intensity[, idx_M0, drop = FALSE]
  I_M1 <- pkm$intensity[, idx_M1, drop = FALSE]

  # Observed ratios
  R_obs_pixel <- I_M1 / pmax(I_M0, 1e-6)
  if (ratio_method == "mean") {
    R_obs <- colMeans(R_obs_pixel, na.rm = TRUE)
  } else if (ratio_method == "median") {
    R_obs <- apply(R_obs_pixel, 2, median, na.rm = TRUE)
  } else {
    R_obs <- colSums(I_M1, na.rm = TRUE) / colSums(I_M0, na.rm = TRUE)
  }

  # Theoretical CIR (M+1 only)
  mz_mono <- pkm$mass[idx_M0]
  R_theo <- cir_ratio(mz_mono)

  # CIR validation
  cir_rel_err <- abs(R_obs - R_theo) / pmax(R_theo, 1e-6)
  is_valid_c13 <- cir_rel_err < cir_rel_tol
  cir_score <- pmin(1, 1 - cir_rel_err)  # [0,1]

  # data frame as result
  data.frame(
    idx_M0 = idx_M0,
    R_obs = R_obs,
    R_theo = R_theo,
    cir_rel_err = cir_rel_err,
    is_valid_c13 = is_valid_c13,
    cir_score = cir_score
  )
}

#apply function: cir_results <-result_cir<- cir_score(result, pkm, min_score_final = 0.6)
