test_that("citrate and isocitrate (same formula) are never fused into one identity", {
  # C6H8O7 is shared by several distinct compounds in the example
  # compound_db, including "Citric acid" (CHEBI 30769) and "Isocitric acid"
  # (CHEBI 30887). Only "Citric acid" shares a ChEBI id with the HCCA
  # standard "CITRATE"; formula alone must never fuse the two compounds.
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
    matrix = "HCCA",
    ppm_tol = 10,
    # Many compounds share the exact mass 192.027; retain all of them so the
    # test isn't sensitive to arbitrary tie-breaking within top_n.
    top_n = NULL,
    quiet = TRUE
  )

  citric <- result[result$candidate_name %in% "Citric acid", ]
  isocitric <- result[result$candidate_name %in% "Isocitric acid", ]

  expect_gt(nrow(citric), 0)
  expect_gt(nrow(isocitric), 0)

  # Distinct rows: formula-only overlap must not merge them.
  expect_true(all(citric$has_standard_compound_match %in% TRUE))
  expect_true(all(isocitric$has_standard_compound_match %in% FALSE))
  expect_true(all(is.na(isocitric$standard_matched_name)))
  expect_true(all(citric$standard_matched_name %in% "CITRATE"))
})


test_that("resolve_compound_identity collapses duplicate compounds across sources", {
  # Mirrors compound_db rows for InChIKey ADDZHRRCUWNSCS-UHFFFAOYSA-N,
  # reported once under Source == "HMDB" and once under Source == "NORMAN".
  compound_db_test <- data.frame(
    Source = c("HMDB", "NORMAN"),
    DB_ID = c("HMDB0001", "NORMAN0001"),
    Name = c("2-Benzofurancarboxaldehyde", "Benzofuran-2-carboxaldehyde"),
    MolecularFormula = c("C9H6O2", "C9H6O2"),
    MonoisotopicMass = c(146.0368, 146.0368),
    StdInChI = NA_character_,
    StdInChIKey = c("ADDZHRRCUWNSCS-UHFFFAOYSA-N", "ADDZHRRCUWNSCS-UHFFFAOYSA-N"),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  resolved <- resolve_compound_identity(compound_db_test)
  expect_equal(length(unique(resolved$compound_identity_id)), 1)

  deduped <- deduplicate_compound_identity(resolved)
  expect_equal(nrow(deduped), 1)
})


test_that("duplicate compound_db rows sharing an InChIKey collapse before matching", {
  compound_db_test <- data.frame(
    Source = c("HMDB", "NORMAN"),
    DB_ID = c("HMDB0001", "NORMAN0001"),
    Name = c("2-Benzofurancarboxaldehyde", "Benzofuran-2-carboxaldehyde"),
    MolecularFormula = c("C9H6O2", "C9H6O2"),
    MonoisotopicMass = c(146.0368, 146.0368),
    StdInChI = NA_character_,
    StdInChIKey = c("ADDZHRRCUWNSCS-UHFFFAOYSA-N", "ADDZHRRCUWNSCS-UHFFFAOYSA-N"),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L, neutral_mass_consensus = 146.0368,
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = 1L, idx = 1L, mz = 147.044076, adduct = "[M+H]+",
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = 1L,
    mz = 147.044076,
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
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 10,
    quiet = TRUE
  )

  candidates <- result[!is.na(result$candidate_name), ]
  expect_equal(nrow(candidates), 1)
})


test_that("compute_standard_adduct_recovery is NA without a standards_db match", {
  expect_true(is.na(compute_standard_adduct_recovery(
    NA_character_, c("[M+H]+"), default_adducts("pos")$name
  )))
})


test_that("compute_standard_adduct_recovery ignores structurally undetectable adduct forms", {
  testable <- default_adducts("pos")$name

  # standards_db lists a dimer and a sodium dimer alongside the one
  # monomer form the pipeline can actually produce - only [M+H]+ should
  # count towards the score, so full recovery of the testable part must
  # read 1, not 1/3.
  score_mixed <- compute_standard_adduct_recovery(
    "[M+H]+;[2M+H]+;[2M+Na]+", c("[M+H]+"), testable
  )
  expect_equal(score_mixed, 1)

  # standards_db lists only forms the pipeline cannot produce at all -
  # this is an impossible comparison, not zero recovery.
  score_untestable_only <- compute_standard_adduct_recovery(
    "[2M+H]+;[2M+Na]+", c("[M+H]+"), testable
  )
  expect_true(is.na(score_untestable_only))

  # every listed adduct is testable - behaviour must match the
  # unfiltered calculation (no regression for the common case).
  score_all_testable <- compute_standard_adduct_recovery(
    "[M+H]+;[M+Na]+", c("[M+H]+"), testable
  )
  expect_equal(score_all_testable, 0.5)
})


test_that("standard_adduct_recovery_score is not deflated by dimer forms end-to-end", {
  data("example_pkm", package = "PeakGuideR")

  res <- run_peakguider_workflow(
    pkm = example_pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE
  )

  baseline <- res$candidate_annotations |>
    dplyr::filter(standard_db_name == "PHENOL", !is.na(standard_adduct_recovery_score))
  expect_true(nrow(baseline) > 0)
  expect_true(all(baseline$standard_adduct_recovery_score == 1))
  expect_equal(baseline$standard_db_adducts[1], "[M+H]+")

  standards_db_dimer <- load_standards_adduct_library(quiet = TRUE)
  phenol_row <- standards_db_dimer[standards_db_dimer$Master_List_NAME == "PHENOL" &
                                      standards_db_dimer$adduct == "[M+H]+", ][1, ]
  expect_equal(nrow(phenol_row), 1)

  dimer_row <- phenol_row
  dimer_row$adduct <- "[2M+H]+"
  standards_db_dimer <- rbind(standards_db_dimer, dimer_row)

  compound_db_default <- load_compound_mass_database(quiet = TRUE)

  res_dimer <- build_candidate_annotations(
    adduct_fam = res$adduct_families,
    feature_summary = res$feature_summary,
    pkm = example_pkm,
    ion_mode = "pos",
    matrix = "HCCA",
    compound_db = compound_db_default,
    standards_db = standards_db_dimer,
    include_single_adduct = TRUE,
    quiet = TRUE
  )

  with_dimer <- res_dimer |>
    dplyr::filter(standard_db_name == "PHENOL", !is.na(standard_adduct_recovery_score))
  expect_true(nrow(with_dimer) > 0)
  listed_adducts <- unique(unlist(strsplit(with_dimer$standard_db_adducts, ";", fixed = TRUE)))
  expect_setequal(listed_adducts, c("[M+H]+", "[2M+H]+"))

  # With the bug, this would read 0.5 (1 of 2 listed forms) instead of 1
  # (1 of 1 testable form) - the phantom dimer must not deflate it.
  expect_true(all(with_dimer$standard_adduct_recovery_score == 1))
})


test_that("standard_adduct_recovery_score/eips_evidence_score = NA are excluded from priority_score, not treated as 0", {
  weights <- default_priority_weights()

  # isotope_evidence_score is left NA on purpose: unlike
  # standard_adduct_recovery_score/eips_evidence_score, it is universally
  # evaluable evidence, so it is zero-filled (kept in the weighted average)
  # rather than excluded - this test isolates the two components that are
  # still excluded/renormalized.
  scores_no_std <- data.frame(
    mass_error_score = 0.9,
    adduct_spatial_score = 0.8,
    isotope_evidence_score = NA_real_,
    eips_evidence_score = NA_real_,
    standard_adduct_recovery_score = NA_real_,
    family_coherence_score = 0.7
  )

  score_na <- compute_priority_score(scores_no_std, weights)

  expected_weight_sum <- weights["mass_error_score"] +
    weights["adduct_spatial_score"] + weights["isotope_evidence_score"] +
    weights["family_coherence_score"]
  expected_score <- (
    0.9 * weights["mass_error_score"] +
      0.8 * weights["adduct_spatial_score"] +
      0 * weights["isotope_evidence_score"] +
      0.7 * weights["family_coherence_score"]
  ) / expected_weight_sum

  expect_equal(unname(score_na), unname(expected_score), tolerance = 1e-8)

  scores_zero <- scores_no_std
  scores_zero$standard_adduct_recovery_score <- 0
  scores_zero$eips_evidence_score <- 0
  score_zero <- compute_priority_score(scores_zero, weights)

  expect_false(isTRUE(all.equal(unname(score_na), unname(score_zero))))
})


test_that("mass_error_score/adduct_spatial_score/isotope_evidence_score/family_coherence_score = NA are zero-filled in priority_score, not excluded", {
  weights <- default_priority_weights()

  scores_mass_only <- data.frame(
    mass_error_score = 0.9,
    adduct_spatial_score = NA_real_,
    isotope_evidence_score = NA_real_,
    eips_evidence_score = NA_real_,
    standard_adduct_recovery_score = NA_real_,
    family_coherence_score = NA_real_
  )

  score_zero_filled <- compute_priority_score(scores_mass_only, weights)

  # Only standard_adduct_recovery_score/eips_evidence_score are excluded and
  # renormalized away; the other three NA components are treated as 0 and
  # keep their weight in the denominator.
  expected_weight_sum <- weights["mass_error_score"] +
    weights["adduct_spatial_score"] + weights["isotope_evidence_score"] +
    weights["family_coherence_score"]
  expected_score <- (0.9 * weights["mass_error_score"]) / expected_weight_sum

  expect_equal(unname(score_zero_filled), unname(expected_score), tolerance = 1e-8)

  # Old behaviour would have excluded every NA component and renormalized
  # down to mass_error_score alone, giving a score equal to mass_error_score
  # itself - i.e. a mass-only candidate could reach ~0.9. The new score must
  # be substantially lower, since three universally-evaluable evidence
  # components are now real negative evidence rather than inapplicable.
  old_style_score <- 0.9
  expect_lt(unname(score_zero_filled), old_style_score - 0.5)
})


test_that("isomeric standards_db compounds within tolerance are scored independently", {
  standards_db_test <- data.frame(
    COMPOUND_ID = c("S1", "S1", "S2"),
    Master_List_NAME = c("Isomer A", "Isomer A", "Isomer B"),
    POLARITY = c("pos", "pos", "pos"),
    matrix = rep("HCCA-DEA solid ionic matrix", 3),
    adduct = c("[M+H]+", "[M+Na]+", "[M+H]+"),
    MOLECULAR_FORMULA = rep("C6H12O6", 3),
    NEUTRAL_MONOISOTOPIC_MASS = rep(180.0634, 3),
    HMDB_clean = NA_character_,
    ChEBI = NA_character_,
    SMILES = NA_character_,
    InCHIKey = c(
      "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
      "BBBBBBBBBBBBBB-UHFFFAOYSA-N"
    ),
    stringsAsFactors = FALSE
  )

  neutral_masses <- data.frame(neutral_mass_id = 1L, neutral_mass_consensus = 180.0634)

  matches <- match_standards_by_mass(
    neutral_masses, standards_db_test,
    ion_mode = "pos", matrix = "HCCA", ppm_tol = 5
  )

  expect_equal(nrow(matches), 2)
  expect_setequal(matches$candidate_name, c("Isomer A", "Isomer B"))

  detected <- c("[M+H]+")
  testable <- default_adducts("pos")$name
  rec_a <- compute_standard_adduct_recovery(
    matches$standard_db_adducts[matches$candidate_name == "Isomer A"], detected, testable
  )
  rec_b <- compute_standard_adduct_recovery(
    matches$standard_db_adducts[matches$candidate_name == "Isomer B"], detected, testable
  )

  expect_equal(rec_a, 0.5)
  expect_equal(rec_b, 1)
  expect_false(isTRUE(all.equal(rec_a, rec_b)))
})


test_that("ambiguous_isomeric flags near-tied top candidates within a neutral_mass_id", {
  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = c(1L, 2L),
      neutral_mass_consensus = c(200.000, 300.000),
      mean_score_adduct = c(0.8, 0.8),
      median_score_adduct = c(0.8, 0.8),
      has_role_conflict = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    family_members = data.frame(
      family_id = c(1L, 2L),
      idx = c(1L, 2L),
      mz = c(201.007276, 301.007276),
      adduct = c("[M+H]+", "[M+H]+"),
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = c(1L, 2L),
    mz = c(201.007276, 301.007276),
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
    DB_ID = c("A1", "A2", "B1"),
    Name = c("Isomer 1", "Isomer 2", "Unique compound"),
    MolecularFormula = c("C10H20O2", "C10H20O2", "C15H30O2"),
    MonoisotopicMass = c(200.0000, 200.0000, 300.0000),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
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

  tied_nmid <- ann$neutral_mass_id[ann$broad_db_name %in% "Isomer 1"]
  tied <- ann[ann$neutral_mass_id %in% tied_nmid, ]
  expect_equal(nrow(tied), 2)
  expect_true(all(tied$ambiguous_isomeric))

  unique_row <- ann[ann$broad_db_name %in% "Unique compound", ]
  expect_equal(nrow(unique_row), 1)
  expect_true(all(!unique_row$ambiguous_isomeric))
})
