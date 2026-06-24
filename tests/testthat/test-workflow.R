test_that("workflow runs on example data", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  expect_s3_class(res, "peakguider_workflow")

  expect_true(all(c(
    "morph_results",
    "cir_results",
    "eips_results",
    "adduct_edges",
    "adduct_families",
    "relation_table",
    "feature_summary",
    "neutral_mass_candidates",
    "pkm",
    "parameters"
  ) %in% names(res)))

  expect_true(is.data.frame(res$feature_summary))
  expect_true(is.data.frame(res$neutral_mass_candidates))
  expect_true(is.data.frame(res$relation_table))
})
