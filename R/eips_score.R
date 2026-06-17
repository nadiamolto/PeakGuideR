#' Score elemental isotope-pattern candidates
#'
#' @description
#' Evaluates non-carbon isotope-pattern candidates by comparing the observed
#' isotope ratio against theoretical element-count ratios for N, O, S, Cl and Br.
#'
#' The function can optionally add formula support by comparing inferred neutral
#' masses with a database-derived EIPS lookup table.
#'
#' Two optional reliability gates can be applied:
#' \itemize{
#'   \item For N/O/S, EIPS is only evaluated if the same `idx_M0` has a
#'         CIR-validated C13_M1 isotope.
#'   \item EIPS is not evaluated for pairs already detected as C13_M2 for an
#'         `idx_M0` with CIR-validated C13_M1.
#' }
#'
#' For low-abundance N/O/S isotope channels, the function can require independent
#' C13 support for the same monoisotopic feature. This reduces spurious elemental
#' isotope assignments in noisy MSI peak tables.
#'
#' In addition, inferred atom counts (`n_hat`) are bounded by the maximum number
#' of atoms observed in the database for that element. Candidates exceeding that
#' bound are discarded before final selection.
#'
#' Expected isotope peak order for the elemental channel evaluated here.
#' N is evaluated at M+1; O, S, Cl and Br are evaluated at M+2.
#'
#' EIPS does not perform the initial mass-difference search. It evaluates
#' elemental isotope candidates already detected by `iso_morphology_candidates()`.
#'
#' @param result Output of `iso_morphology_candidates()` (must contain at least
#'   `idx_M0`, `idx_cand`, `score_final`, `element`, `iso_type`).
#' @param pkm Peak matrix list with `mass` and `intensity`.
#' @param eips_n_table Optional precomputed theoretical lookup table with
#'   columns `element`, `k`, `delta`, `n` and `R_theo`. If `NULL`, the internal
#'   PeakGuideR `eips_n_table` object is used.
#' @param eips_table Optional database-derived lookup table with columns
#'   `mz_mono`, `element` and `n_el`. If `NULL`, the internal PeakGuideR
#'   `eips_table` object is used.
#' @param ion_mode `"pos"` or `"neg"` (used only when formula support is enabled).
#' @param adducts Optional data.frame with adduct definitions. If `NULL`,
#'   `default_adducts()` is used. The table must contain `name`, `mode` and
#'   `mass`, where `mass` is the net m/z shift relative to the neutral mass.
#' @param min_score_final Minimum morphology score required to consider a pair.
#' @param eips_rel_tol Relative tolerance used as validity threshold and for
#'   ratio score scaling.
#' @param ratio_method `"sum"`, `"mean"`, or `"median"` aggregation of isotope ratios.
#'   `"sum"` computes sum(I_iso) / sum(I_M0) over selected pixels and is recommended
#'   for MSI isotope-ratio estimation.
#' @param min_quantile Quantile used for intensity masking. Default= 0.01.
#' @param ppm_neutral Ppm tolerance used when evaluating neutral mass
#' consistency for elemental isotope-pattern support.
#' @param n_window Allowed deviation between inferred `n_hat` and formula `n_el`
#'   when formula support is evaluated.
#' @param top_n_hat Integer. Keep top-k best `n_hat` candidates by relative error.
#' @param cir_df Optional CIR results table. Must contain `idx_M0` and
#'   `is_valid_c13`. If supplied, enables the N/O/S gate.
#' @param morph_df Optional morphology table (typically the same object passed in
#'   `result`). Must contain `idx_M0`, `idx_cand`, `iso_type`. If supplied,
#'   enables exclusion of C13_M2 pairs.
#' @param require_c13_for Character vector of elements for which CIR-validated
#'   C13_M1 is required. Default: `c("N","O","S")`.
#' @param exclude_c13_m2 Logical. If `TRUE`, excludes pairs already detected as
#'   C13_M2 for `idx_M0` values with CIR-valid C13_M1.
#' @param return_debug Logical. If `TRUE`, returns a list with final output and
#'   internal top-n candidates.
#'
#' @return By default, a single `data.frame` / tibble with one row per evaluated pair:
#' \itemize{
#'   \item `idx_M0`, `mz_M0`, `idx_cand`, `mz_iso`
#'   \item `iso_type`, `element`, `k`, `delta`
#'   \item `R_obs`, `n_hat`, `R_theo_hat`, `eips_rel_err`
#'   \item `is_valid_eips`, `score_eips`
#'   \item `has_formula_support`, `adduct_name`, `neutral_mass`,
#'         `mz_mono_db`, `n_el_db`, `neutral_err_ppm`
#' }
#'
#' If `return_debug = TRUE`, returns a list with:
#' \itemize{
#'   \item `eips_validation`
#'   \item `top_n_candidates`
#' }
#'
#' @export
eips_score <- function(result, pkm,
                       eips_n_table=NULL,
                       eips_table = NULL,
                       ion_mode = c("pos","neg"),
                       adducts = NULL,
                       min_score_final = 0.6,
                       eips_rel_tol = 0.3,
                       ratio_method = c("sum","median","mean"),
                       min_quantile = 0.01,
                       ppm_neutral = 5,
                       n_window = 1L,
                       top_n_hat = 5L,
                       cir_df = NULL,
                       morph_df = NULL,
                       require_c13_for = c("N","O","S"),
                       exclude_c13_m2 = TRUE,
                       return_debug = FALSE) {

  ion_mode     <- match.arg(ion_mode)
  ratio_method <- match.arg(ratio_method)

  if (is.null(eips_n_table)) {
    eips_n_table <- get("eips_n_table", envir = asNamespace("PeakGuideR"))
  }

  gaussian_score <- function(err, tol) {
    err <- pmax(0, err)
    exp(-(err^2) / (2 * tol^2))
  }

  if (is.null(adducts)) adducts <- default_adducts()

  required_result_cols <- c("idx_M0", "idx_cand", "score_final", "element", "iso_type")

  missing_result_cols <- setdiff(required_result_cols, names(result))
  if (length(missing_result_cols) > 0) {
    stop(
      "`result` is missing required columns: ",
      paste(missing_result_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.list(pkm) || is.null(pkm$mass) || is.null(pkm$intensity)) {
    stop("`pkm` must be a list containing `mass` and `intensity`.", call. = FALSE)
  }

  if (!is.numeric(pkm$mass)) {
    stop("`pkm$mass` must be a numeric vector.", call. = FALSE)
  }

  if (!is.matrix(pkm$intensity) || !is.numeric(pkm$intensity)) {
    stop("`pkm$intensity` must be a numeric matrix.", call. = FALSE)
  }

  if (length(pkm$mass) != ncol(pkm$intensity)) {
    stop(
      "Length of `pkm$mass` must match the number of columns in `pkm$intensity`.",
      call. = FALSE
    )
  }

  stopifnot(
    is.data.frame(eips_n_table),
    all(c("element", "k", "delta", "n", "R_theo") %in% names(eips_n_table))
  )

  stopifnot(
    is.data.frame(adducts),
    all(c("name", "mode", "mass") %in% names(adducts))
  )
  if (is.null(eips_table)) {
    eips_table <- get("eips_table", envir = asNamespace("PeakGuideR"))
  }

  stopifnot(is.data.frame(eips_table),
            all(c("mz_mono","element","n_el") %in% names(eips_table)))

  top_n_hat <- as.integer(top_n_hat)
  if (!is.finite(top_n_hat) || top_n_hat < 1L) top_n_hat <- 5L


  # Hard upper bound per element from DB
  max_n_db <- eips_table |>
    dplyr::group_by(element) |>
    dplyr::summarise(n_max_db = max(n_el, na.rm = TRUE), .groups = "drop")

  # Clean theoretical table
  eips_n_theo <- eips_n_table |>
    dplyr::transmute(
      element,
      k_theo      = as.integer(k),
      delta_theo  = delta,
      n_theo      = as.integer(n),
      R_theo_theo = as.numeric(R_theo)
    )

  # 1. Filter EIPS candidates from morphology
  pairs <- result |>
    dplyr::filter(!is.na(score_final), score_final >= min_score_final) |>
    dplyr::filter(!is.na(element), element %in% c("N","O","S","Cl","Br"))


  if (nrow(pairs) == 0) {
    empty <- dplyr::tibble()
    return(if (isTRUE(return_debug)) list(eips_validation = empty,
                                          top_n_candidates = empty) else empty)
  }

  if (any(pairs$idx_M0 < 1 | pairs$idx_M0 > ncol(pkm$intensity)) ||
      any(pairs$idx_cand < 1 | pairs$idx_cand > ncol(pkm$intensity))) {
    stop("Feature indices in `result` are outside the columns of `pkm$intensity`.", call. = FALSE)
  }
  # 2. Reliability gates based on CIR-valid C13_M1
  valid_c13_m0 <- integer(0)
  if (!is.null(cir_df) && is.data.frame(cir_df) &&
      all(c("idx_M0","is_valid_c13") %in% names(cir_df))) {
    valid_c13_m0 <- cir_df |>
      dplyr::filter(is_valid_c13 %in% TRUE) |>
      dplyr::pull(idx_M0) |>
      unique()
  }

  # Gate A: N/O/S only evaluated if the same M0 has a CIR-validated C13_M1
  pairs <- pairs |>
    dplyr::filter(
      !(element %in% require_c13_for) |
        idx_M0 %in% valid_c13_m0
    )

  # Gate B: exclude pairs already labelled as C13_M2 under CIR-valid M0
  if (isTRUE(exclude_c13_m2) &&
      length(valid_c13_m0) > 0 &&
      !is.null(morph_df) && is.data.frame(morph_df) &&
      all(c("idx_M0","idx_cand","iso_type") %in% names(morph_df))) {

    c13_m2_reserved <- morph_df |>
      dplyr::filter(idx_M0 %in% valid_c13_m0, iso_type == "C13_M2") |>
      dplyr::distinct(idx_M0, idx_cand)

    pairs <- pairs |>
      dplyr::anti_join(c13_m2_reserved, by = c("idx_M0","idx_cand"))
  }

  if (nrow(pairs) == 0) {
    empty <- dplyr::tibble()
    return(if (isTRUE(return_debug)) list(eips_validation = empty,
                                          top_n_candidates = empty) else empty)
  }

  # 3. Observed ratios
  idx_M0  <- pairs$idx_M0
  idx_iso <- pairs$idx_cand

  I_M0  <- pkm$intensity[, idx_M0,  drop = FALSE]
  I_iso <- pkm$intensity[, idx_iso, drop = FALSE]

  min_I_M0  <- apply(I_M0,  2, function(x) stats::quantile(x, probs = min_quantile, na.rm = TRUE))
  min_I_iso <- apply(I_iso, 2, function(x) stats::quantile(x, probs = min_quantile, na.rm = TRUE))

  # Element-dependent masking
  mask_strategy <- dplyr::case_when(
    pairs$element %in% c("N","O","S") ~ "M0_only",
    pairs$element %in% c("Cl","Br")   ~ "AND",
    TRUE                              ~ "AND"
  )

  mask_valid <- matrix(FALSE, nrow = nrow(I_M0), ncol = ncol(I_M0))
  for (j in seq_len(ncol(I_M0))) {
    if (mask_strategy[j] == "M0_only") {
      mask_valid[, j] <- I_M0[, j] > min_I_M0[j]
    } else {
      mask_valid[, j] <- (I_M0[, j] > min_I_M0[j]) & (I_iso[, j] > min_I_iso[j])
    }
  }

  I_M0_f  <- I_M0
  I_iso_f <- I_iso
  I_M0_f[!mask_valid]  <- NA
  I_iso_f[!mask_valid] <- NA

  R_obs_pixel <- I_iso_f / pmax(I_M0_f, 1e-6)

  R_obs <- if (ratio_method == "sum") {
    colSums(I_iso_f, na.rm = TRUE) / pmax(colSums(I_M0_f, na.rm = TRUE), 1e-6)
  } else if (ratio_method == "mean") {
    colMeans(R_obs_pixel, na.rm = TRUE)
  } else {
    apply(R_obs_pixel, 2, stats::median, na.rm = TRUE)
  }

  pairs2 <- pairs |>
    dplyr::mutate(
      mz_M0  = pkm$mass[idx_M0],
      mz_iso = pkm$mass[idx_iso],
      R_obs  = as.numeric(R_obs)
    )

  # Expected channel per element
  element_k <- c(N = 1L, O = 2L, S = 2L, Cl = 2L, Br = 2L)

  pairs2 <- pairs2 |>
    dplyr::mutate(k_expected = unname(element_k[element]))

  # 4. Ratio-only theoretical matching
  joined_n <- pairs2 |>
    dplyr::inner_join(eips_n_theo, by = "element", relationship = "many-to-many") |>
    dplyr::filter(k_theo == k_expected) |>
    dplyr::mutate(
      eips_rel_err = abs(R_obs - R_theo_theo) / pmax(R_theo_theo, 1e-6)
    ) |>
    dplyr::left_join(max_n_db, by = "element") |>
    dplyr::mutate(
      n_max_db = dplyr::if_else(is.na(n_max_db), Inf, n_max_db),
      is_valid_eips_bound = n_theo <= n_max_db
    ) |>
    dplyr::filter(is_valid_eips_bound)   # <- if above DB max, it does not appear

  if (nrow(joined_n) == 0) {
    empty <- dplyr::tibble()
    return(if (isTRUE(return_debug)) list(eips_validation = empty,
                                          top_n_candidates = empty) else empty)
  }

  top_n_candidates <- joined_n |>
    dplyr::group_by(idx_M0, idx_cand, element) |>
    dplyr::slice_min(order_by = eips_rel_err, n = top_n_hat, with_ties = TRUE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_hat         = n_theo,
      R_theo_hat    = R_theo_theo,
      is_valid_eips = eips_rel_err < eips_rel_tol,
      score_eips    = gaussian_score(eips_rel_err, eips_rel_tol)
    )

  # Keep only best ratio candidate per pair
  best_ratio <- top_n_candidates |>
    dplyr::group_by(idx_M0, idx_cand, element) |>
    dplyr::slice_min(order_by = eips_rel_err, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx_M0, mz_M0,
      idx_cand, mz_iso,
      iso_type, element,
      k = k_expected,
      delta = dplyr::coalesce(delta_theo, NA_real_),
      R_obs,
      n_hat,
      R_theo_hat,
      eips_rel_err,
      is_valid_eips,
      score_eips
    )

  # 5. Optional formula support

  adf <- adducts |>
    dplyr::filter(mode == ion_mode)

  if (!nrow(adf)) {
    stop("No adducts available for ion_mode = '", ion_mode, "'.", call. = FALSE)
  }
  expd <- tidyr::crossing(
    best_ratio,
    adf |> dplyr::select(adduct_name = name, adduct_mass = mass)
  ) |>
    dplyr::mutate(
      neutral_mass = mz_M0 - adduct_mass
    )

  supp <- expd |>
    dplyr::inner_join(eips_table, by = "element", relationship = "many-to-many") |>
    dplyr::filter(
      abs(neutral_mass - mz_mono) <= (mz_mono * ppm_neutral * 1e-6),
      abs(as.integer(n_el) - as.integer(n_hat)) <= as.integer(n_window)
    ) |>
    dplyr::mutate(
      neutral_err_ppm = 1e6 * abs(neutral_mass - mz_mono) / pmax(mz_mono, 1e-6)
    ) |>
    dplyr::group_by(idx_M0, idx_cand, element) |>
    dplyr::slice_min(order_by = neutral_err_ppm, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      idx_M0, idx_cand, element,
      has_formula_support = TRUE,
      adduct_name,
      neutral_mass,
      mz_mono_db = mz_mono,
      n_el_db = as.integer(n_el),
      neutral_err_ppm
    )


  # 6. Single final table
  eips_validation <- best_ratio |>
    dplyr::left_join(supp, by = c("idx_M0","idx_cand","element")) |>
    dplyr::mutate(
      has_formula_support = dplyr::coalesce(has_formula_support, FALSE)
    )

  if (isTRUE(return_debug)) {
    return(list(
      eips_validation = eips_validation,
      top_n_candidates = top_n_candidates
    ))
  } else {
    return(eips_validation)
  }
}


#eips_results <- eips_score( result = res_morph, pkm = pkm, eips_n_table = eips_n_table, eips_table = eips_table,
#ion_mode = "pos", min_score_final = 0.4, eips_rel_tol = 0.3, ppm_neutral = 10, n_window = 2L,
#ratio_method = "sum",cir_df = cir_results, morph_df = res_morph)
