test_that("broad_db_only confidence_class reflects the number of supporting signals", {
  compound_db_test <- data.frame(
    Source = rep("test_db", 3),
    DB_ID = c("C1", "C2", "C3"),
    Name = c("Mass Only Compound", "Single Evidence Compound", "Multi Evidence Compound"),
    MolecularFormula = c("C6H10O2", "C7H12O2", "C8H14O2"),
    MonoisotopicMass = c(114.0681, 128.0837, 142.0994),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = c(1L, 2L, 3L),
      neutral_mass_consensus = c(114.0681, 128.0837, 142.0994),
      mean_score_adduct = c(0.1, 0.1, 0.1),
      median_score_adduct = c(0.1, 0.1, 0.1),
      has_role_conflict = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = c(1L, 2L, 3L),
      idx = c(1L, 2L, 3L),
      mz = c(115.075376, 129.091018, 143.106658),
      adduct = c("[M+H]+", "[M+H]+", "[M+H]+"),
      stringsAsFactors = FALSE
    )
  )

  # adduct_spatial_score = 0.1 in every family (below the default
  # adduct_min_score = 0.5), so only is_c13_m0/has_eips drive the count:
  # 0 signals, exactly 1, and both (2) signals.
  feature_summary_test <- data.frame(
    idx = c(1L, 2L, 3L),
    mz = c(115.075376, 129.091018, 143.106658),
    is_c13_m0 = c(FALSE, TRUE, TRUE),
    c13_score = c(NA_real_, 0.9, 0.9),
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = c(FALSE, FALSE, TRUE),
    eips_elements = NA_character_,
    eips_score = c(NA_real_, NA_real_, 0.9),
    stringsAsFactors = FALSE
  )

  ann <- build_candidate_annotations(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    pkm = NULL,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 10,
    top_n = NULL,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  by_name <- stats::setNames(ann$confidence_class, ann$broad_db_name)

  expect_true(all(ann$source %in% "broad_db_only"))
  expect_equal(unname(by_name["Mass Only Compound"]), "broad_db_only_mass_only")
  expect_equal(unname(by_name["Single Evidence Compound"]), "broad_db_only_single_evidence")
  expect_equal(unname(by_name["Multi Evidence Compound"]), "broad_db_only_multi_evidence")
})


test_that("cross-database-linked confidence_class is split by standard-adduct recovery", {
  # "Test Standard A"/"Test Standard B" each annotate 3 adducts in
  # standards_db - restricted to names actually present in
  # default_adducts("neg") ([M-H]-, [M+Cl]-, [M-H-H2O]-), since adducts
  # standards_db lists that the detection pipeline could never produce
  # (e.g. the "[M+Na-2H]-"/"[M+K-2H]-" this fixture used before) are
  # excluded from standard_adduct_recovery_score's denominator and would
  # not exercise the intended 3/3-vs-1/3 contrast below. The family for A
  # recovers 3/3 (>= default recovery_threshold = 0.5) and the family for
  # B recovers 1/3 (< 0.5).
  compound_db_test <- data.frame(
    Source = c("test_db", "test_db"),
    DB_ID = c("A1", "B1"),
    Name = c("Test Compound A", "Test Compound B"),
    MolecularFormula = c("C10H10O2", "C15H20O3"),
    MonoisotopicMass = c(200.0000, 300.0000),
    StdInChI = NA_character_,
    StdInChIKey = c("AAAAAAAAAAAAAA-UHFFFAOYSA-N", "BBBBBBBBBBBBBB-UHFFFAOYSA-N"),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  standards_db_test <- data.frame(
    COMPOUND_ID = c(rep("S1", 3), rep("S2", 3)),
    Master_List_NAME = c(rep("Test Standard A", 3), rep("Test Standard B", 3)),
    POLARITY = "neg",
    matrix = "HCCA-DEA solid ionic matrix",
    adduct = rep(c("[M-H]-", "[M+Cl]-", "[M-H-H2O]-"), 2),
    MOLECULAR_FORMULA = c(rep("C10H10O2", 3), rep("C15H20O3", 3)),
    NEUTRAL_MONOISOTOPIC_MASS = c(rep(200.0000, 3), rep(300.0000, 3)),
    HMDB_clean = NA_character_,
    ChEBI = NA_character_,
    SMILES = NA_character_,
    InCHIKey = c(
      rep("AAAAAAAAAAAAAA-UHFFFAOYSA-N", 3),
      rep("BBBBBBBBBBBBBB-UHFFFAOYSA-N", 3)
    ),
    stringsAsFactors = FALSE
  )

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = c(1L, 2L),
      neutral_mass_consensus = c(200.0000, 300.0000),
      mean_score_adduct = c(0.9, 0.9),
      median_score_adduct = c(0.9, 0.9),
      has_role_conflict = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = c(1L, 1L, 1L, 2L),
      idx = c(1L, 2L, 3L, 4L),
      mz = c(201.007276, 222.989218, 238.963158, 301.007276),
      adduct = c("[M-H]-", "[M+Cl]-", "[M-H-H2O]-", "[M-H]-"),
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = c(1L, 2L, 3L, 4L),
    mz = c(201.007276, 222.989218, 238.963158, 301.007276),
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
    compound_db = compound_db_test,
    standards_db = standards_db_test,
    ion_mode = "neg",
    matrix = "HCCA",
    ppm_tol = 10,
    top_n = NULL,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  by_name <- stats::setNames(ann$confidence_class, ann$broad_db_name)

  expect_true(all(ann$source %in% "both"))
  expect_equal(unname(by_name["Test Compound A"]), "cross_db_linked_adduct_recovered")
  expect_equal(unname(by_name["Test Compound B"]), "cross_db_linked_partial_recovery")
})
