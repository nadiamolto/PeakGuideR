test_that("eips_score returns a data frame", {
  data("example_pkm", package = "PeakGuideR")

  morph <- iso_morphology_candidates(
    pkm = example_pkm,
    prefer_mode = "ppm",
    method = "pearson",
    transform = "none"
  )

  cir <- cir_score(
    result = morph,
    pkm = example_pkm,
    min_score_final = 0.6,
    cir_rel_tol = 0.3,
    ratio_method = "sum"
  )

  eips <- eips_score(
    result = morph,
    pkm = example_pkm,
    ion_mode = "pos",
    adducts = default_adducts("pos"),
    min_score_final = 0.6,
    eips_rel_tol = 0.3,
    ratio_method = "sum",
    cir_df = cir,
    morph_df = morph
  )

  expect_true(is.data.frame(eips))
})
