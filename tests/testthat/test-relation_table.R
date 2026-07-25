
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

test_that("build_relation_table() returns the full schema when there is no evidence at all", {
  rt <- build_relation_table(cir_results = NULL, eips_results = NULL, adduct_fam = NULL)

  expected_cols <- c(
    "from_idx", "to_idx", "from_mz", "to_mz",
    "relation_type", "evidence_type", "evidence_score", "score_type",
    "is_valid", "group_id", "from_role", "to_role", "label"
  )

  expect_equal(nrow(rt), 0L)
  expect_true(all(expected_cols %in% names(rt)))
  expect_type(rt$from_idx, "integer")
  expect_type(rt$to_idx, "integer")
  expect_type(rt$from_mz, "double")
  expect_type(rt$to_mz, "double")
  expect_type(rt$relation_type, "character")
  expect_type(rt$evidence_type, "character")
  expect_type(rt$evidence_score, "double")
  expect_type(rt$score_type, "character")
  expect_type(rt$is_valid, "logical")
  expect_type(rt$group_id, "integer")
  expect_type(rt$from_role, "character")
  expect_type(rt$to_role, "character")
  expect_type(rt$label, "character")

  # Same empty-evidence case with a data.frame that has 0 rows (as opposed
  # to fully NULL inputs), to also cover the "!length(edges)" branch when
  # cir_results/eips_results are empty data.frames rather than NULL.
  empty_cir <- data.frame(
    idx_M0 = integer(0), idx_M1 = integer(0),
    mz_M0 = numeric(0), mz_M1 = numeric(0),
    cir_score = numeric(0), is_valid_c13 = logical(0)
  )
  rt2 <- build_relation_table(cir_results = empty_cir, eips_results = NULL, adduct_fam = NULL)
  expect_equal(nrow(rt2), 0L)
  expect_true(all(expected_cols %in% names(rt2)))

  # build_feature_summary() must run without error on the empty result -
  # this is exactly what broke before the fix.
  fs <- build_feature_summary(rt)
  expect_s3_class(fs, "data.frame")
  expect_equal(nrow(fs), 0L)
})
