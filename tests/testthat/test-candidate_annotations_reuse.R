test_that("build_neutral_mass_candidates() gives the same result with a pre-computed candidate_annotations", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  old_result <- build_neutral_mass_candidates(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  precomputed_annotations <- build_candidate_annotations(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    include_single_adduct = TRUE,
    quiet = TRUE
  )

  new_result <- build_neutral_mass_candidates(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    ion_mode = "pos",
    matrix = "HCCA",
    candidate_annotations = precomputed_annotations,
    quiet = TRUE
  )

  key_cols <- c("neutral_mass_id", "candidate_ppm_error", "candidate_name", "candidate_db_id")

  old_sorted <- old_result |> dplyr::arrange(dplyr::across(dplyr::all_of(key_cols)))
  new_sorted <- new_result |> dplyr::arrange(dplyr::across(dplyr::all_of(key_cols)))

  expect_equal(old_sorted, new_sorted)
})


test_that("build_neutral_mass_candidates() only keeps hypothesis_origin == \"family\" from a supplied candidate_annotations", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  # Some neutral_mass_id values are shared between a family and a
  # single_adduct_clustered hypothesis; without filtering to
  # hypothesis_origin == "family" first, those rows would be matched twice
  # against the same compound/standard candidates and duplicated.
  expect_true(any(
    res$candidate_annotations$hypothesis_origin %in% "single_adduct_clustered"
  ))

  result <- build_neutral_mass_candidates(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    ion_mode = "pos",
    matrix = "HCCA",
    candidate_annotations = res$candidate_annotations,
    quiet = TRUE
  )

  reference <- build_neutral_mass_candidates(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  expect_equal(nrow(result), nrow(reference))

  dup_key <- paste(
    result$neutral_mass_id, result$candidate_source,
    result$candidate_db_id, result$candidate_name,
    sep = "__"
  )
  expect_equal(anyDuplicated(dup_key), 0L)
})


test_that("run_peakguider_workflow(include_single_adduct = FALSE) excludes single-adduct hypotheses", {
  data("example_pkm", package = "PeakGuideR")

  res_no_single <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  expect_false(any(
    res_no_single$candidate_annotations$hypothesis_origin %in%
      c("single_adduct_isolated", "single_adduct_clustered")
  ))

  res_default <- run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  )

  expect_true(any(
    res_default$candidate_annotations$hypothesis_origin %in%
      c("single_adduct_isolated", "single_adduct_clustered")
  ))
})


test_that("run_peakguider_workflow() computes candidate_annotations only once", {
  data("example_pkm", package = "PeakGuideR")

  ns <- asNamespace("PeakGuideR")

  options(peakguider_test_build_candidate_annotations_calls = 0L)
  on.exit(options(peakguider_test_build_candidate_annotations_calls = NULL), add = TRUE)

  suppressMessages(
    trace(
      "build_candidate_annotations",
      tracer = quote({
        options(
          peakguider_test_build_candidate_annotations_calls =
            getOption("peakguider_test_build_candidate_annotations_calls", 0L) + 1L
        )
      }),
      print = FALSE,
      where = ns
    )
  )
  on.exit(
    suppressMessages(
      try(untrace("build_candidate_annotations", where = ns), silent = TRUE)
    ),
    add = TRUE
  )

  invisible(run_peakguider_workflow(
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    quiet = TRUE
  ))

  expect_equal(getOption("peakguider_test_build_candidate_annotations_calls"), 1L)
})
