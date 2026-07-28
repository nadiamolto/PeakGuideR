# Minimal hand-built pkm shared by the tests below: a 3x3 pixel grid,
# 3 features (an M0/M+1 isotope-like pair at idx 1/2, plus an unrelated
# adduct-like feature at idx 3).
make_test_pkm <- function() {
  pos <- cbind(x = rep(1:3, 3), y = rep(1:3, each = 3))
  intensity <- matrix(
    c(
      10, 20, 30, 40, 50, 60, 70, 80, 90,   # feature 1 (M0)
      1,  2,  3,  4,  5,  6,  7,  8,  9,    # feature 2 (M+1)
      100, 90, 80, 70, 60, 50, 40, 30, 20   # feature 3
    ),
    nrow = 9, ncol = 3
  )
  list(mass = c(500.0000, 501.0034, 522.9814), intensity = intensity, pos = pos)
}

test_that("plot_ion_image returns a ggplot object with a valid idx", {
  skip_if_not_installed("ggplot2")
  pkm <- make_test_pkm()

  p <- plot_ion_image(pkm, idx = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ion_image errors when both idx and mz are supplied", {
  skip_if_not_installed("ggplot2")
  pkm <- make_test_pkm()

  expect_error(plot_ion_image(pkm, idx = 1, mz = 500), "Only one of")
})

test_that("plot_ion_image errors when neither idx nor mz is supplied", {
  skip_if_not_installed("ggplot2")
  pkm <- make_test_pkm()

  expect_error(plot_ion_image(pkm), "Exactly one of")
})

test_that("plot_ion_image finds the closest feature by mz within tolerance", {
  skip_if_not_installed("ggplot2")
  pkm <- make_test_pkm()

  p <- plot_ion_image(pkm, mz = 500.0005, tol_ppm = 10)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, sprintf("m/z %.4f", pkm$mass[1]))
})

test_that("plot_ion_image errors when no feature is within tol_ppm of mz", {
  skip_if_not_installed("ggplot2")
  pkm <- make_test_pkm()

  expect_error(plot_ion_image(pkm, mz = 600, tol_ppm = 10), "No feature found")
})

test_that("plot_isotope_pair returns a patchwork object and adds an M+2 panel when available", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pkm <- make_test_pkm()

  morph_row <- data.frame(
    idx_M0 = 1L, mz_M0 = pkm$mass[1], idx_cand = 2L, mz_cand = pkm$mass[2],
    score_final = 0.87, iso_type = "C13_M1", element = "C"
  )
  p_morph <- plot_isotope_pair(pkm, idx_M0 = 1L, idx_cand = 2L, morph_row = morph_row)
  expect_s3_class(p_morph, "patchwork")

  eips_row <- data.frame(
    idx_M0 = 1L, mz_M0 = pkm$mass[1], idx_cand = 3L, mz_iso = pkm$mass[3],
    element = "S", R_obs = 0.05, n_hat = 1, R_theo_hat = 0.045,
    score_eips = 0.9, is_valid_eips = TRUE
  )
  p_eips <- plot_isotope_pair(pkm, idx_M0 = 1L, idx_cand = 3L, eips_row = eips_row)
  expect_s3_class(p_eips, "patchwork")

  cir_row_no_m2 <- data.frame(
    idx_M0 = 1L, mz_M0 = pkm$mass[1], idx_M1 = 2L, mz_M1 = pkm$mass[2],
    R_obs = 0.108, R_theo = 0.110, cir_score = 0.95, cir_class = "high_agreement",
    is_valid_c13 = TRUE, has_C13_M2 = FALSE, idx_C13_M2 = NA_integer_,
    mz_C13_M2 = NA_real_, score_C13_M2 = NA_real_
  )
  p_cir <- plot_isotope_pair(pkm, idx_M0 = 1L, idx_cand = 2L, cir_row = cir_row_no_m2)
  expect_s3_class(p_cir, "patchwork")
  expect_equal(length(p_cir$patches$plots), 1)  # base plot + 1 extra = 2 panels total

  cir_row_with_m2 <- cir_row_no_m2
  cir_row_with_m2$has_C13_M2 <- TRUE
  cir_row_with_m2$idx_C13_M2 <- 3L
  cir_row_with_m2$mz_C13_M2 <- pkm$mass[3]
  cir_row_with_m2$score_C13_M2 <- 0.8

  p_cir_m2 <- plot_isotope_pair(pkm, idx_M0 = 1L, idx_cand = 2L, cir_row = cir_row_with_m2)
  expect_s3_class(p_cir_m2, "patchwork")
  expect_equal(length(p_cir_m2$patches$plots), 2)  # base plot + 2 extra = 3 panels total
})

test_that("plot_isotope_pair validates row-argument columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pkm <- make_test_pkm()

  bad_morph_row <- data.frame(idx_M0 = 1L, mz_M0 = pkm$mass[1])
  expect_error(
    plot_isotope_pair(pkm, idx_M0 = 1L, idx_cand = 2L, morph_row = bad_morph_row),
    "missing required columns"
  )
})

test_that("plot_adduct_family returns a patchwork object for a hand-built family", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pkm <- make_test_pkm()

  adduct_fam <- list(
    family_summary = data.frame(family_id = 1L, neutral_mass_consensus = 499.0),
    family_members = data.frame(
      family_id = c(1L, 1L),
      idx = c(1L, 2L),
      mz = pkm$mass[c(1, 2)],
      adduct = c("[M+H]+", "[M+Na]+"),
      neutral_mass_from_feature = c(499.0, 499.0),
      n_edges = c(1L, 1L),
      mean_edge_score = c(0.9, 0.9)
    )
  )

  p <- plot_adduct_family(pkm, family_id = 1L, adduct_fam = adduct_fam)
  expect_s3_class(p, "patchwork")
  expect_equal(length(p$patches$plots), 1)  # base plot + 1 extra = 2 members total
})

test_that("plot_adduct_family errors on an unknown family_id", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pkm <- make_test_pkm()

  adduct_fam <- list(
    family_summary = data.frame(family_id = 1L, neutral_mass_consensus = 499.0),
    family_members = data.frame(
      family_id = 1L, idx = 1L, mz = pkm$mass[1], adduct = "[M+H]+",
      neutral_mass_from_feature = 499.0, n_edges = 1L, mean_edge_score = 0.9
    )
  )

  expect_error(
    plot_adduct_family(pkm, family_id = 99L, adduct_fam = adduct_fam),
    "not found"
  )
})
