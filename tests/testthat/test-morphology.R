test_that("iso_morphology_candidates returns expected columns", {
  data("example_pkm", package = "PeakGuideR")

  morph <- iso_morphology_candidates(
    pkm = example_pkm,
    prefer_mode = "ppm",
    method = "pearson",
    transform = "none"
  )

  expect_true(is.data.frame(morph))

  expected_cols <- c(
    "idx_M0",
    "mz_M0",
    "iso_type",
    "element",
    "k",
    "z",
    "idx_cand",
    "mz_cand",
    "score_global",
    "score_final",
    "mass_err_da",
    "mass_err_ppm",
    "mass_dev_score"
  )

  expect_true(all(expected_cols %in% names(morph)))
})
