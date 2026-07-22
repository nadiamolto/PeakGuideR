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


test_that("top_n limits candidates per neutral mass", {

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = 100,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L,
      idx = 1L,
      mz = 101.0073,
      adduct = "[M+H]+",
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 101.0073,
    is_c13_m0 = FALSE,
    c13_score = NA_real_,
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = FALSE,
    eips_elements = NA_character_,
    eips_score = NA_real_,
    stringsAsFactors = FALSE
  )

  compound_db_test <- data.frame(
    Source = rep("test_db", 3),
    DB_ID = paste0("ID", 1:3),
    Name = paste0("Compound_", 1:3),
    MolecularFormula = c("C1", "C2", "C3"),
    MonoisotopicMass = c(
      100.0000,
      100.0001,
      100.0002
    ),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- build_neutral_mass_candidates(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 5,
    top_n = 2,
    quiet = TRUE
  )

  candidates <- result |>
    dplyr::filter(!is.na(candidate_name))

  expect_equal(nrow(candidates), 2)
  expect_equal(
    candidates$candidate_name,
    c("Compound_1", "Compound_2")
  )
})

test_that("top_n NULL retains all candidates within tolerance", {

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = 100,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L,
      idx = 1L,
      mz = 101.0073,
      adduct = "[M+H]+",
      stringsAsFactors = FALSE
    )
  )
  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 101.0073,
    is_c13_m0 = FALSE,
    c13_score = NA_real_,
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = FALSE,
    eips_elements = NA_character_,
    eips_score = NA_real_,
    stringsAsFactors = FALSE
  )


  compound_db_test <- data.frame(
    Source = rep("test_db", 3),
    DB_ID = paste0("ID", 1:3),
    Name = paste0("Compound_", 1:3),
    MolecularFormula = c("C1", "C2", "C3"),
    MonoisotopicMass = c(
      100.0000,
      100.0001,
      100.0002
    ),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- build_neutral_mass_candidates(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 5,
    top_n = NULL,
    quiet = TRUE
  )

  candidates <- result |>
    dplyr::filter(!is.na(candidate_name))

  expect_equal(nrow(candidates), 3)
  expect_equal(
    candidates$candidate_name,
    c("Compound_1", "Compound_2", "Compound_3")
  )

  expect_true(
    all(candidates$candidate_ppm_error <= 5)
  )
})


test_that("top_n must be NULL or a positive integer", {

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = 100,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L,
      idx = 1L,
      mz = 101.0073,
      adduct = "[M+H]+",
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 101.0073,
    stringsAsFactors = FALSE
  )

  compound_db_test <- data.frame(
    Source = "test_db",
    DB_ID = "ID1",
    Name = "Compound_1",
    MolecularFormula = "C1",
    MonoisotopicMass = 100,
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  expect_error(
    build_neutral_mass_candidates(
      adduct_fam = adduct_fam_test,
      feature_summary = feature_summary_test,
      compound_db = compound_db_test,
      standards_db = NULL,
      ion_mode = "pos",
      matrix = NULL,
      top_n = 0,
      quiet = TRUE
    ),
    "`top_n` must be NULL or a positive integer",
    fixed = TRUE
  )
})

