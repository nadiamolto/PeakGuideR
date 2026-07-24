test_that("match_standards_by_mass() with matrix = NULL skips standards_db entirely", {
  standards_db <- load_standards_adduct_library(quiet = TRUE)
  neutral_masses <- data.frame(
    neutral_mass_id = 1L,
    neutral_mass_consensus = standards_db$NEUTRAL_MONOISOTOPIC_MASS[1]
  )

  expect_message(
    res <- match_standards_by_mass(
      neutral_masses, standards_db,
      ion_mode = "neg", matrix = NULL
    ),
    "matrix"
  )
  expect_equal(nrow(res), 0)

  expect_silent(
    match_standards_by_mass(
      neutral_masses, standards_db,
      ion_mode = "neg", matrix = NULL, quiet = TRUE
    )
  )
})


test_that("match_standards_by_mass() with matrix = NULL never touches standards_db", {
  neutral_masses <- data.frame(neutral_mass_id = 1L, neutral_mass_consensus = 192.027)

  res <- match_standards_by_mass(
    neutral_masses, standards_db = NULL,
    ion_mode = "neg", matrix = NULL, quiet = TRUE
  )

  expect_equal(nrow(res), 0)
})


test_that("match_standards_by_mass() with an unrelated matrix name returns no matches", {
  standards_db <- load_standards_adduct_library(quiet = TRUE)
  neutral_masses <- data.frame(
    neutral_mass_id = 1L,
    neutral_mass_consensus = standards_db$NEUTRAL_MONOISOTOPIC_MASS[1]
  )

  res <- match_standards_by_mass(
    neutral_masses, standards_db,
    ion_mode = "neg", matrix = "DHB", quiet = TRUE
  )

  expect_equal(nrow(res), 0)
})


test_that("match_standards_by_mass() with matrix = \"HCCA\" is unchanged", {
  standards_db <- load_standards_adduct_library(quiet = TRUE)
  neutral_masses <- data.frame(
    neutral_mass_id = 1L,
    neutral_mass_consensus = standards_db$NEUTRAL_MONOISOTOPIC_MASS[1]
  )

  res <- match_standards_by_mass(
    neutral_masses, standards_db,
    ion_mode = "neg", matrix = "HCCA", quiet = TRUE
  )

  expect_equal(nrow(res), 1)
  expect_equal(res$candidate_name, "CITRATE")
  expect_equal(res$candidate_source, "standards_db")
})


test_that("build_candidate_annotations() with matrix = NULL never produces standards-linked candidates", {
  # C6H8O7 (citric acid neutral mass) also matches the "CITRATE" standard by
  # mass when matrix = "HCCA" is used; with matrix = NULL, standards_db must
  # not be consulted at all, regardless of what would otherwise match.
  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = 192.02700,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L,
      idx = 1L,
      mz = 191.019724,
      adduct = "[M-H]-",
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 191.019724,
    is_c13_m0 = FALSE,
    c13_score = NA_real_,
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = FALSE,
    eips_elements = NA_character_,
    eips_score = NA_real_,
    stringsAsFactors = FALSE
  )

  ann <- build_candidate_annotations(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    pkm = NULL,
    ion_mode = "neg",
    matrix = NULL,
    ppm_tol = 10,
    top_n = NULL,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  expect_gt(nrow(ann), 0)
  expect_true(all(ann$source %in% "broad_db_only"))
  expect_false(any(ann$source %in% c("both", "standards_only")))
})


test_that("build_neutral_mass_candidates() with matrix = NULL never sets has_standard_compound_match", {
  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = 192.02700,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L,
      idx = 1L,
      mz = 191.019724,
      adduct = "[M-H]-",
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 191.019724,
    is_c13_m0 = FALSE,
    c13_score = NA_real_,
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = FALSE,
    eips_elements = NA_character_,
    eips_score = NA_real_,
    stringsAsFactors = FALSE
  )

  result <- build_neutral_mass_candidates(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    ion_mode = "neg",
    matrix = NULL,
    ppm_tol = 10,
    top_n = NULL,
    quiet = TRUE
  )

  expect_gt(nrow(result), 0)
  expect_false(any(result$has_standard_compound_match %in% TRUE))
  expect_false(any(result$has_standard_adduct_match %in% TRUE))
})
