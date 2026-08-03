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


single_adduct_isolated_fixture <- function(compound_db_test, mz, intensity) {
  adduct_fam_empty <- list(
    family_summary = data.frame(
      family_id = integer(0), neutral_mass_consensus = numeric(0),
      mean_score_adduct = numeric(0), median_score_adduct = numeric(0),
      has_role_conflict = logical(0)
    ),
    family_members = data.frame(
      family_id = integer(0), idx = integer(0), mz = numeric(0), adduct = character(0),
      stringsAsFactors = FALSE
    )
  )

  feature_summary_test <- data.frame(
    idx = seq_along(mz),
    mz = mz,
    is_c13_m0 = FALSE, c13_score = NA_real_, has_c13_m2_support = FALSE, c13_m2_score = NA_real_,
    has_eips = FALSE, eips_elements = NA_character_, eips_score = NA_real_,
    has_adduct_family = FALSE, is_c13_m1 = FALSE, is_c13_m2 = FALSE, is_eips_isotope = FALSE,
    stringsAsFactors = FALSE
  )

  pkm_test <- list(mass = mz, intensity = intensity)

  build_candidate_annotations(
    adduct_fam = adduct_fam_empty,
    feature_summary = feature_summary_test,
    pkm = pkm_test,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 10,
    neutral_cluster_ppm = 10,
    top_n = NULL,
    include_single_adduct = TRUE,
    quiet = TRUE
  )
}


test_that("ambiguous_isomeric is not falsely triggered by two single_adduct_isolated rows of the same identity", {
  # Reproduces the Toufik brain dataset pattern: [M+Na]+ and [M+K]+ of the
  # same compound are each inferred independently, without ever forming a
  # shared adduct family. Both hypotheses cluster into the same
  # neutral_mass_id and both resolve to the same compound_db candidate, so
  # they carry an identical priority_score - a naive row-vs-row comparison
  # reads that as a tied top-1/top-2 and flags ambiguity where there is only
  # one identity.
  compound_db_test <- data.frame(
    Source = "test_db",
    DB_ID = "C1",
    Name = "Test Compound",
    MolecularFormula = "C10H20O2",
    MonoisotopicMass = 200.000000,
    StdInChI = NA_character_,
    StdInChIKey = "AAAAAAAAAAAAAA-BBBBBBBBBB-A",
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  mz <- c(200.000000 + 22.989218, 200.000000 + 38.963158) # [M+Na]+, [M+K]+
  intensity <- matrix(
    c(10, 12, 9, 11, 13, 8, 14, 10, 9, 12,
      11, 13, 8, 12, 14, 9, 13, 11, 8, 13),
    nrow = 10, ncol = 2
  )

  ann <- single_adduct_isolated_fixture(compound_db_test, mz, intensity)

  dup <- ann[!is.na(ann$broad_db_name) & ann$broad_db_name == "Test Compound", ]
  expect_equal(nrow(dup), 2)
  expect_setequal(dup$hypothesis_origin, "single_adduct_isolated")
  expect_equal(length(unique(dup$neutral_mass_id)), 1)
  expect_false(any(dup$ambiguous_isomeric))
})


test_that("ambiguous_isomeric still flags genuinely competing identities under single_adduct_isolated", {
  # Same two-adduct, no-family setup as above, but with two distinct
  # compound_db entries at the same mass instead of one: each identity is
  # still duplicated across its own [M+Na]+/[M+K]+ rows, and that
  # duplication must not mask the real ambiguity between Identity A and
  # Identity B.
  compound_db_test <- data.frame(
    Source = c("test_db", "test_db"),
    DB_ID = c("A1", "B1"),
    Name = c("Identity A", "Identity B"),
    MolecularFormula = c("C10H20O2", "C10H20O2"),
    MonoisotopicMass = c(200.000000, 200.000000),
    StdInChI = NA_character_,
    StdInChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-A", "CCCCCCCCCCCCCC-DDDDDDDDDD-C"),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  mz <- c(200.000000 + 22.989218, 200.000000 + 38.963158) # [M+Na]+, [M+K]+
  intensity <- matrix(
    c(10, 12, 9, 11, 13, 8, 14, 10, 9, 12,
      11, 13, 8, 12, 14, 9, 13, 11, 8, 13),
    nrow = 10, ncol = 2
  )

  ann <- single_adduct_isolated_fixture(compound_db_test, mz, intensity)

  sub <- ann[!is.na(ann$broad_db_name), ]
  expect_equal(nrow(sub), 4)
  expect_equal(length(unique(sub$neutral_mass_id)), 1)
  expect_true(all(sub$ambiguous_isomeric))
})


test_that("ambiguous_isomeric compares a 3-row identity against the next distinct identity, not its own rows", {
  # Identity A is recovered via three independent single-adduct hypotheses
  # ([M+H]+, [M+Na]+, [M+K]+) that never merged into a family, so it
  # occupies three rows of its own for this neutral_mass_id. Identity B
  # occupies another three rows at a slightly worse mass error. Under the
  # old row-vs-row comparison, Identity A's own top two rows are tied
  # (identical priority_score), which alone was enough to flag the mass as
  # ambiguous regardless of Identity B. After collapsing by identity, the
  # comparison is A's best row vs B's best row, with a gap well above
  # ambiguity_gap - not ambiguous.
  compound_db_test <- data.frame(
    Source = c("test_db", "test_db"),
    DB_ID = c("A1", "B1"),
    Name = c("Identity A", "Identity B"),
    MolecularFormula = c("C10H20O2", "C10H20O2"),
    MonoisotopicMass = c(200.000000, 200.001600),
    StdInChI = NA_character_,
    StdInChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-A", "CCCCCCCCCCCCCC-DDDDDDDDDD-C"),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  mz <- c(
    200.000000 + 1.007276, # [M+H]+
    200.000000 + 22.989218, # [M+Na]+
    200.000000 + 38.963158 # [M+K]+
  )
  intensity <- matrix(
    c(10, 12, 9, 11, 13, 8, 14, 10, 9, 12,
      11, 13, 8, 12, 14, 9, 13, 11, 8, 13,
      9, 11, 10, 13, 12, 8, 11, 10, 9, 14),
    nrow = 10, ncol = 3
  )

  ann <- single_adduct_isolated_fixture(compound_db_test, mz, intensity)

  identity_a <- ann[!is.na(ann$broad_db_name) & ann$broad_db_name == "Identity A", ]
  identity_b <- ann[!is.na(ann$broad_db_name) & ann$broad_db_name == "Identity B", ]
  expect_equal(nrow(identity_a), 3)
  expect_equal(nrow(identity_b), 3)
  expect_equal(length(unique(c(identity_a$neutral_mass_id, identity_b$neutral_mass_id))), 1)

  # sanity check: Identity A's own rows are tied at the top, so a
  # row-vs-row comparison (the old behaviour) would flag this ambiguous
  # regardless of Identity B.
  expect_equal(sum(identity_a$priority_score == max(identity_a$priority_score)), 2)

  expect_false(any(identity_a$ambiguous_isomeric))
  expect_false(any(identity_b$ambiguous_isomeric))
})


test_that("expand_compound_candidates_for_standards() recovers a top_n-truncated standard-linked candidate", {
  # 5 compound_db isobars crowd the same mass window; "Isobar 4" is the one
  # that shares an InChIKey with a standards_db compound, but it has the
  # worst (4th-ranked) ppm error of the five, so a plain top_n = 3 cut
  # truncates it away before resolve_cross_source_identity() ever sees it.
  compound_db_test <- data.frame(
    Source = rep("test_db", 5),
    DB_ID = c("C1", "C2", "C3", "C4", "C5"),
    Name = c("Isobar 1", "Isobar 2", "Isobar 3", "Isobar 4", "Isobar 5"),
    MolecularFormula = rep("C5H10O2", 5),
    MonoisotopicMass = c(100.000000, 100.000100, 100.000200, 100.000300, 100.000400),
    StdInChI = NA_character_,
    StdInChIKey = c(NA, NA, NA, "STDSTDSTDSTDST-UHFFFAOYSA-S", NA),
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  standards_db_test <- data.frame(
    COMPOUND_ID = "S1",
    Master_List_NAME = "Isobar 4 standard",
    POLARITY = "pos",
    matrix = "HCCA-DEA solid ionic matrix",
    adduct = "[M+H]+",
    MOLECULAR_FORMULA = "C5H10O2",
    NEUTRAL_MONOISOTOPIC_MASS = 100.000300,
    HMDB_clean = NA_character_,
    ChEBI = NA_character_,
    SMILES = NA_character_,
    InCHIKey = "STDSTDSTDSTDST-UHFFFAOYSA-S",
    stringsAsFactors = FALSE
  )

  neutral_masses <- data.frame(neutral_mass_id = 1L, neutral_mass_consensus = 100.000000)

  compound_candidates <- match_compound_db_by_mass(
    neutral_masses, compound_db_test, ppm_tol = 10, top_n = 3
  )
  expect_equal(nrow(compound_candidates), 3)
  expect_false("Isobar 4" %in% compound_candidates$candidate_name)

  standard_candidates <- match_standards_by_mass(
    neutral_masses, standards_db_test, ion_mode = "pos", matrix = "HCCA", ppm_tol = 10
  )
  expect_equal(nrow(standard_candidates), 1)

  expanded <- expand_compound_candidates_for_standards(
    compound_candidates, standard_candidates, neutral_masses, compound_db_test, 10
  )
  expect_true("Isobar 4" %in% expanded$candidate_name)

  resolved <- resolve_cross_source_identity(expanded, standard_candidates)
  isobar4_row <- resolved[resolved$broad_db_name %in% "Isobar 4", ]
  expect_equal(nrow(isobar4_row), 1)
  expect_equal(isobar4_row$source, "both")
  expect_equal(isobar4_row$standard_db_name, "Isobar 4 standard")
})


test_that("expand_compound_candidates_for_standards() leaves top_n untouched when no standard matches", {
  compound_db_test <- data.frame(
    Source = rep("test_db", 5),
    DB_ID = c("C1", "C2", "C3", "C4", "C5"),
    Name = c("Isobar 1", "Isobar 2", "Isobar 3", "Isobar 4", "Isobar 5"),
    MolecularFormula = rep("C5H10O2", 5),
    MonoisotopicMass = c(100.000000, 100.000100, 100.000200, 100.000300, 100.000400),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  neutral_masses <- data.frame(neutral_mass_id = 1L, neutral_mass_consensus = 100.000000)

  compound_candidates <- match_compound_db_by_mass(
    neutral_masses, compound_db_test, ppm_tol = 10, top_n = 3
  )
  expect_equal(nrow(compound_candidates), 3)

  expanded_no_std <- expand_compound_candidates_for_standards(
    compound_candidates, NULL, neutral_masses, compound_db_test, 10
  )
  expect_equal(expanded_no_std, compound_candidates)

  standard_candidates_empty <- compound_candidates[0, ]
  expanded_empty_std <- expand_compound_candidates_for_standards(
    compound_candidates, standard_candidates_empty, neutral_masses, compound_db_test, 10
  )
  expect_equal(expanded_empty_std, compound_candidates)
})
