test_that("internal reference tables are available without attaching the package", {
  ns <- asNamespace("PeakGuideR")

  expect_true(exists("cir_table", envir = ns, inherits = FALSE))
  expect_true(exists("eips_table", envir = ns, inherits = FALSE))
  expect_true(exists("eips_n_table", envir = ns, inherits = FALSE))

  exports <- getNamespaceExports("PeakGuideR")

  expect_false("cir_table" %in% exports)
  expect_false("eips_table" %in% exports)
  expect_false("eips_n_table" %in% exports)
})
