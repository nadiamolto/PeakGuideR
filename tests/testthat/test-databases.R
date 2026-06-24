test_that("example databases can be loaded", {
  compound_db_path <- system.file(
    "extdata/compound_mass_database_example.rds",
    package = "PeakGuideR"
  )

  standards_db_path <- system.file(
    "extdata/standards_adduct_library_example.rds",
    package = "PeakGuideR"
  )

  expect_true(file.exists(compound_db_path))
  expect_true(file.exists(standards_db_path))

  compound_db <- readRDS(compound_db_path)
  standards_db <- readRDS(standards_db_path)

  expect_true(is.data.frame(compound_db))
  expect_true(is.data.frame(standards_db))
  expect_gt(nrow(compound_db), 0)
  expect_gt(nrow(standards_db), 0)
})
