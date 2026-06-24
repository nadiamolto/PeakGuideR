test_that("neutral mass candidates table has expected columns", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  required_cols <- c(
    "neutral_mass_id",
    "neutral_mass_consensus",
    "candidate_name",
    "candidate_ppm_error"
  )

  expect_true(all(required_cols %in% names(res$neutral_mass_candidates)))
})
