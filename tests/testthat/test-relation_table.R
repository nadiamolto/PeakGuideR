
test_that("relation table contains required columns", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  required_cols <- c(
    "from_idx",
    "to_idx",
    "from_mz",
    "to_mz",
    "relation_type",
    "evidence_type",
    "evidence_score",
    "from_role",
    "to_role"
  )

  expect_true(all(required_cols %in% names(res$relation_table)))
})
