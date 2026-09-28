test_that("build_neutral_mass_candidates() keeps feature_idx/feature_mz/inferred_adducts positionally aligned", {
  # Regression test: feature_idx, feature_mz and inferred_adducts used to be
  # built with two independent sort() calls (one numeric on idx/mz, one
  # alphabetical on the adduct label), so position i in inferred_adducts did
  # not necessarily correspond to position i in feature_idx/feature_mz unless
  # the alphabetical order of the labels happened to match the mass order.
  # These 6 adducts (the same set used to diagnose the bug against real data)
  # sort very differently by label than by m/z, so this reproduces it
  # reliably. Rows are deliberately supplied out of idx/mz order too, so the
  # fix's arrange() - not incidental input order - is what's being tested.
  adduct_masses <- c(
    "[M+H-H2O]+" = -17.003289,
    "[M+H]+"     = 1.007276,
    "[M+Na]+"    = 22.989218,
    "[M+K]+"     = 38.963158,
    "[M+2Na-H]+" = 44.971160,
    "[M+2K-H]+"  = 76.919040
  )
  neutral_mass <- 189.0426
  shuffle <- c(4, 1, 6, 2, 5, 3)

  family_members_test <- data.frame(
    family_id = 1L,
    idx = (1:6)[shuffle],
    mz = round(neutral_mass + unname(adduct_masses)[shuffle], 6),
    adduct = names(adduct_masses)[shuffle],
    stringsAsFactors = FALSE
  )

  adduct_fam_test <- list(
    family_summary = data.frame(
      family_id = 1L,
      neutral_mass_consensus = neutral_mass,
      stringsAsFactors = FALSE
    ),
    family_members = family_members_test
  )

  feature_summary_test <- data.frame(
    idx = 1:6,
    mz = round(neutral_mass + unname(adduct_masses), 6),
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
    Source = "test_db",
    DB_ID = "ID1",
    Name = "Compound_1",
    MolecularFormula = "C1",
    MonoisotopicMass = neutral_mass,
    StdInChI = NA_character_,
    StdInChIKey = NA_character_,
    SMILES = NA_character_,
    Kegg = NA_character_,
    stringsAsFactors = FALSE
  )

  check_positional_correspondence <- function(feature_idx_str, feature_mz_str, adduct_str) {
    idx_vals <- as.integer(strsplit(feature_idx_str, ";")[[1]])
    mz_vals <- as.numeric(strsplit(feature_mz_str, ";")[[1]])
    adduct_vals <- strsplit(adduct_str, ";")[[1]]

    expect_equal(length(idx_vals), 6L)
    expect_equal(length(mz_vals), 6L)
    expect_equal(length(adduct_vals), 6L)

    # idx should come out ascending (mirrors ascending m/z here)
    expect_equal(idx_vals, sort(idx_vals))

    implied_neutral <- mz_vals - unname(adduct_masses[adduct_vals])
    ppm_error <- 1e6 * abs(implied_neutral - neutral_mass) / neutral_mass
    expect_true(all(ppm_error < 0.6))
  }

  result <- build_neutral_mass_candidates(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    ppm_tol = 5,
    top_n = 1,
    quiet = TRUE
  )

  row <- result[1, ]
  check_positional_correspondence(row$feature_idx, row$feature_mz, row$inferred_adducts)

  annotations <- build_candidate_annotations(
    adduct_fam = adduct_fam_test,
    feature_summary = feature_summary_test,
    compound_db = compound_db_test,
    standards_db = NULL,
    ion_mode = "pos",
    matrix = NULL,
    top_n = 1,
    include_single_adduct = FALSE,
    quiet = TRUE
  )

  fam_row <- annotations[annotations$hypothesis_origin == "family", ][1, ]
  feature_mz_from_idx <- feature_summary_test$mz[
    as.integer(strsplit(fam_row$feature_idx, ";")[[1]])
  ]
  check_positional_correspondence(
    fam_row$feature_idx,
    paste(feature_mz_from_idx, collapse = ";"),
    fam_row$inferred_adduct
  )
})
