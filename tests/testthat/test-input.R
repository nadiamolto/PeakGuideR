test_that("workflow fails clearly without pkm", {
  expect_error(
    run_peakguider_workflow(
      ion_mode = "pos",
      matrix = "HCCA",
      quiet = TRUE
    ),
    "pkm"
  )
})

test_that("workflow fails with malformed pkm", {
  bad_pkm <- list(mass = c(100, 101))

  expect_error(
    run_peakguider_workflow(
      pkm = bad_pkm,
      ion_mode = "pos",
      matrix = "HCCA",
      quiet = TRUE
    ),
    "mass.*intensity"
  )
})
