# EIPS evidence that is fully explained by a feature's own assigned adduct
# (a 35Cl/37Cl isotope pattern on a feature annotated `[M+Cl]-`) provides no
# independent support that the neutral molecule contains that element, and must
# not inflate the per-neutral-mass EIPS aggregates. Evidence for elements the
# adduct does not introduce (N/O/S, or Cl on any other adduct) must be scored
# exactly as before.

adduct_masses_neg <- c(
  "[M-H]-"  = -1.007276,
  "[M+Cl]-" = 34.968853
)

make_eips_fixture <- function() {
  # Three families, one neutral mass each:
  #  - 200: Cl-EIPS sits on the `[M+Cl]-` feature  -> adduct-explained, dropped
  #  - 300: N-EIPS sits on the `[M-H]-` feature     -> genuine, untouched
  #  - 400: Cl-EIPS sits on the `[M-H]-` feature    -> genuine, untouched
  neutral <- c(fam_cl = 200.0, fam_n = 300.0, fam_cl_real = 400.0)

  family_members <- data.frame(
    family_id = c(1L, 1L, 2L, 2L, 3L, 3L),
    idx = 1:6,
    adduct = c("[M-H]-", "[M+Cl]-", "[M-H]-", "[M+Cl]-", "[M-H]-", "[M+Cl]-"),
    stringsAsFactors = FALSE
  )
  family_members$mz <- round(
    rep(unname(neutral), each = 2) +
      unname(adduct_masses_neg[family_members$adduct]),
    6
  )

  adduct_fam <- list(
    family_summary = data.frame(
      family_id = 1:3,
      neutral_mass_consensus = unname(neutral),
      stringsAsFactors = FALSE
    ),
    family_members = family_members
  )

  feature_summary <- data.frame(
    idx = 1:6,
    mz = family_members$mz[order(family_members$idx)],
    is_c13_m0 = FALSE,
    c13_score = NA_real_,
    has_c13_m2_support = FALSE,
    c13_m2_score = NA_real_,
    has_eips = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
    eips_elements = c(NA, "Cl", "N", NA, "Cl", NA),
    eips_score = c(NA, 0.90, 0.80, NA, 0.70, NA),
    stringsAsFactors = FALSE
  )

  compound_db <- data.frame(
    Source = "test_db",
    DB_ID = paste0("ID", 1:3),
    Name = paste0("Compound_", 1:3),
    MolecularFormula = c("C1", "C1N1", "C1Cl1"),
    MonoisotopicMass = unname(neutral),
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  list(
    adduct_fam = adduct_fam,
    feature_summary = feature_summary,
    compound_db = compound_db,
    neutral = neutral
  )
}

annotation_for_mass <- function(annotations, neutral_mass, tol = 1e-3) {
  row <- annotations[abs(annotations$neutral_mass_consensus - neutral_mass) < tol, ]
  row[which.max(row$priority_score), ]
}

neutral_row_for_mass <- function(neutral_tbl, neutral_mass, tol = 1e-3) {
  row <- neutral_tbl[abs(neutral_tbl$neutral_mass_consensus - neutral_mass) < tol, ]
  row[1, ]
}


test_that("Cl-EIPS on a `[M+Cl]-` feature does not count as EIPS evidence", {
  fx <- make_eips_fixture()

  annotations <- build_candidate_annotations(
    adduct_fam = fx$adduct_fam,
    feature_summary = fx$feature_summary,
    compound_db = fx$compound_db,
    standards_db = NULL,
    ion_mode = "neg",
    matrix = NULL,
    top_n = 1,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  cl_row <- annotation_for_mass(annotations, fx$neutral[["fam_cl"]])

  expect_true(is.na(cl_row$eips_evidence_score))
  expect_identical(cl_row$confidence_class, "broad_db_only_mass_only")

  neutral_tbl <- build_neutral_mass_candidates(
    adduct_fam = fx$adduct_fam,
    feature_summary = fx$feature_summary,
    compound_db = fx$compound_db,
    standards_db = NULL,
    ion_mode = "neg",
    matrix = NULL,
    top_n = 1,
    quiet = TRUE
  )

  cl_neutral <- neutral_row_for_mass(neutral_tbl, fx$neutral[["fam_cl"]])
  expect_false(cl_neutral$has_EIPS_support)
  expect_true(is.na(cl_neutral$EIPS_score))
  expect_true(is.na(cl_neutral$EIPS_elements))
})


test_that("N-EIPS and Cl-EIPS on non-`[M+Cl]-` features are scored as before", {
  fx <- make_eips_fixture()

  annotations <- build_candidate_annotations(
    adduct_fam = fx$adduct_fam,
    feature_summary = fx$feature_summary,
    compound_db = fx$compound_db,
    standards_db = NULL,
    ion_mode = "neg",
    matrix = NULL,
    top_n = 1,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  n_row <- annotation_for_mass(annotations, fx$neutral[["fam_n"]])
  expect_equal(n_row$eips_evidence_score, 0.80)
  expect_identical(n_row$confidence_class, "broad_db_only_single_evidence")

  cl_real_row <- annotation_for_mass(annotations, fx$neutral[["fam_cl_real"]])
  expect_equal(cl_real_row$eips_evidence_score, 0.70)
  expect_identical(cl_real_row$confidence_class, "broad_db_only_single_evidence")

  neutral_tbl <- build_neutral_mass_candidates(
    adduct_fam = fx$adduct_fam,
    feature_summary = fx$feature_summary,
    compound_db = fx$compound_db,
    standards_db = NULL,
    ion_mode = "neg",
    matrix = NULL,
    top_n = 1,
    quiet = TRUE
  )

  n_neutral <- neutral_row_for_mass(neutral_tbl, fx$neutral[["fam_n"]])
  expect_true(n_neutral$has_EIPS_support)
  expect_equal(n_neutral$EIPS_score, 0.80)
  expect_identical(n_neutral$EIPS_elements, "N")

  cl_real_neutral <- neutral_row_for_mass(neutral_tbl, fx$neutral[["fam_cl_real"]])
  expect_true(cl_real_neutral$has_EIPS_support)
  expect_equal(cl_real_neutral$EIPS_score, 0.70)
  expect_identical(cl_real_neutral$EIPS_elements, "Cl")
})


test_that("drop_adduct_explained_eips() only removes the adduct's own element", {
  # Cl stripped on `[M+Cl]-`; N kept; Cl kept on `[M-H]-`; multi-element
  # feature keeps the element the adduct does not supply.
  out <- drop_adduct_explained_eips(
    adduct = c("[M+Cl]-", "[M-H]-", "[M+Cl]-", "[M+Cl]-"),
    has_eips = c(TRUE, TRUE, TRUE, FALSE),
    eips_elements = c("Cl", "Cl", "Cl;N", NA),
    eips_score = c(0.9, 0.7, 0.8, NA)
  )

  expect_equal(out$has_eips, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(out$eips_elements, c(NA, "Cl", "N", NA))
  expect_equal(out$eips_score, c(NA, 0.7, 0.8, NA))
})
