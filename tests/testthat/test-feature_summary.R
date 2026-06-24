test_that("feature summary has one row per input feature", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  expect_equal(nrow(res$feature_summary), length(example_pkm$mass))
  expect_true(all(c("idx", "mz", "has_adduct_family") %in% names(res$feature_summary)))
})
