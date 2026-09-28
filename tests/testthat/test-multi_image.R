make_two_image_pkm <- function() {
  data("example_pkm", package = "PeakGuideR", envir = environment())
  pkm1 <- example_pkm
  n <- nrow(pkm1$intensity)

  # Second image: same features, different pixel content, and the *same*
  # local pixel-coordinate grid as the first one (coordinates overlap
  # between images, as with separately acquired runs).
  pkm2_intensity <- pkm1$intensity[rev(seq_len(n)), , drop = FALSE]

  pkm <- pkm1
  pkm$intensity <- rbind(pkm1$intensity, pkm2_intensity)
  pkm$pos <- rbind(pkm1$pos, pkm1$pos)
  pkm$posMotors <- rbind(pkm1$posMotors, pkm1$posMotors)
  pkm$SNR <- rbind(pkm1$SNR, pkm1$SNR)
  pkm$area <- rbind(pkm1$area, pkm1$area)
  pkm$numPixels <- c(n, n)
  pkm$names <- c("image_A", "image_B")
  pkm
}


test_that("derive_image_id_from_pkm() expands numPixels/names into a per-pixel factor", {
  pkm <- list(
    intensity = matrix(0, nrow = 5, ncol = 2),
    numPixels = c(2L, 3L),
    names = c("A", "B")
  )

  ids <- derive_image_id_from_pkm(pkm)

  expect_s3_class(ids, "factor")
  expect_equal(as.character(ids), c("A", "A", "B", "B", "B"))
  expect_equal(levels(ids), c("A", "B"))
})


test_that("derive_image_id_from_pkm() prefers an existing image_id", {
  pkm <- list(
    intensity = matrix(0, nrow = 4, ncol = 2),
    numPixels = c(1L, 3L),
    names = c("A", "B"),
    image_id = factor(c("r1", "r1", "r2", "r2"))
  )

  expect_equal(as.character(derive_image_id_from_pkm(pkm)), c("r1", "r1", "r2", "r2"))
})


test_that("derive_image_id_from_pkm() treats single-image input as one image", {
  one <- list(
    intensity = matrix(0, nrow = 4, ncol = 2),
    numPixels = 4L,
    names = "only"
  )
  expect_equal(levels(derive_image_id_from_pkm(one)), "only")
  expect_length(derive_image_id_from_pkm(one), 4L)

  bare <- list(intensity = matrix(0, nrow = 3, ncol = 2))
  expect_equal(nlevels(derive_image_id_from_pkm(bare)), 1L)
  expect_length(derive_image_id_from_pkm(bare), 3L)
})


test_that("derive_image_id_from_pkm() rejects inconsistent inputs", {
  expect_error(
    derive_image_id_from_pkm(list(
      intensity = matrix(0, nrow = 5, ncol = 2),
      numPixels = c(2L, 2L),
      names = c("A", "B")
    )),
    "numPixels"
  )
  expect_error(
    derive_image_id_from_pkm(list(
      intensity = matrix(0, nrow = 5, ncol = 2),
      image_id = c("A", "B")
    )),
    "image_id"
  )
})


test_that("cardinal_to_peakmatrix() captures the per-pixel run as image_id", {
  skip_if_not_installed("Cardinal")
  skip_if_not_installed("BiocParallel")

  suppressPackageStartupMessages(requireNamespace("Cardinal"))
  Cardinal::setCardinalBPPARAM(BiocParallel::SerialParam())

  set.seed(1)
  sim <- Cardinal::simulateImage(
    preset = 1, npeaks = 30, dim = c(6, 5), representation = "centroid"
  )

  # single run
  pkm_single <- cardinal_to_peakmatrix(sim)
  expect_s3_class(pkm_single$image_id, "factor")
  expect_length(pkm_single$image_id, nrow(pkm_single$intensity))
  expect_equal(nlevels(pkm_single$image_id), 1L)

  # two runs
  Cardinal::run(sim) <- factor(rep(c("runA", "runB"), each = ncol(sim) / 2))
  pkm_multi <- cardinal_to_peakmatrix(sim)

  expect_s3_class(pkm_multi$image_id, "factor")
  expect_equal(levels(pkm_multi$image_id), c("runA", "runB"))
  expect_equal(as.integer(table(pkm_multi$image_id)), c(15L, 15L))

  # existing fields are unchanged by the addition
  expect_equal(pkm_multi$numPixels, nrow(pkm_multi$intensity))
  expect_equal(pkm_multi$names, "runA")
  expect_true(all(
    c("mass", "binSize", "intensity", "SNR", "area", "normalizations",
      "pos", "numPixels", "names", "posMotors") %in% names(pkm_multi)
  ))
})


test_that("run_peakguider_workflow(multi_image = TRUE) pools two overlapping images and disables tiles", {
  pkm <- make_two_image_pkm()

  # overlapping coordinates between the two images
  n <- pkm$numPixels[1]
  expect_equal(pkm$pos[seq_len(n), ], pkm$pos[n + seq_len(n), ])

  expect_message(
    res_multi <- run_peakguider_workflow(
      pkm = pkm, ion_mode = "pos", matrix = "HCCA",
      multi_image = TRUE
    ),
    "tile-based isotope-morphology consistency is disabled"
  )

  # (a) runs to completion with one joint result
  expect_s3_class(res_multi, "peakguider_workflow")
  expect_true(is.data.frame(res_multi$candidate_annotations))
  expect_true(is.data.frame(res_multi$neutral_mass_candidates))

  # (b) tile-based step is disabled: no tile statistic, score = global score
  expect_gt(nrow(res_multi$morph_results), 0)
  expect_true(all(is.na(res_multi$morph_results$tile_consistency)))
  expect_true(all(is.na(res_multi$morph_results$tile_summary)))
  expect_equal(
    res_multi$morph_results$score_final,
    res_multi$morph_results$score_global
  )

  # with tiles on (default), the same input does compute tile statistics
  res_tiles <- run_peakguider_workflow(
    pkm = pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE
  )
  expect_true(any(is.finite(res_tiles$morph_results$tile_consistency)))

  # (c) one joint table over all pixels, not split per image
  expect_false("image_id" %in% names(res_multi$candidate_annotations))
  expect_false("image_id" %in% names(res_multi$neutral_mass_candidates))
  expect_false("image_id" %in% names(res_multi$feature_summary))
  expect_equal(anyDuplicated(res_multi$candidate_annotations[
    , c("neutral_mass_id", "feature_idx", "broad_db_id", "standard_db_compound_id",
        "hypothesis_origin")
  ]), 0L)

  # image identity is attached to the returned pkm for plotting
  expect_equal(nrow(res_multi$pkm$intensity), 2 * n)
  expect_equal(levels(res_multi$pkm$image_id), c("image_A", "image_B"))
  expect_true(res_multi$parameters$multi_image)
})


test_that("run_peakguider_workflow() default (multi_image = FALSE) is unchanged", {
  data("example_pkm", package = "PeakGuideR")

  res_default <- run_peakguider_workflow(
    pkm = example_pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE
  )
  res_explicit <- run_peakguider_workflow(
    pkm = example_pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE,
    multi_image = FALSE
  )

  expect_identical(res_default, res_explicit)
  expect_false("multi_image" %in% names(res_default$parameters))
  expect_false("image_id" %in% names(res_default$pkm))
  expect_true(any(is.finite(res_default$morph_results$tile_consistency)))
})


test_that("plot_ion_image() draws a single panel by default and facets by image on request", {
  skip_if_not_installed("ggplot2")
  data("example_pkm", package = "PeakGuideR")
  two <- make_two_image_pkm()

  # single-image input: one panel, no facet, no extra data column
  p_single <- plot_ion_image(example_pkm, idx = 1)
  expect_s3_class(p_single$facet, "FacetNull")
  expect_false("image" %in% names(p_single$data))
  expect_equal(nrow(p_single$data), nrow(example_pkm$intensity))

  # two-image input, image identity auto-derived from numPixels/names
  p_auto <- plot_ion_image(two, idx = 1)
  expect_s3_class(p_auto$facet, "FacetWrap")
  expect_equal(nlevels(p_auto$data$image), 2L)
  expect_equal(levels(p_auto$data$image), c("image_A", "image_B"))
  built <- ggplot2::ggplot_build(p_auto)
  expect_equal(length(unique(built$layout$layout$PANEL)), 2L)

  # explicit grouping vector
  grp <- rep(c("left", "right"), length.out = nrow(two$intensity))
  p_group <- plot_ion_image(two, idx = 1, group = grp)
  expect_equal(nlevels(p_group$data$image), 2L)

  # group = FALSE forces the original single panel
  p_off <- plot_ion_image(two, idx = 1, group = FALSE)
  expect_s3_class(p_off$facet, "FacetNull")

  # a mismatched grouping vector is rejected
  expect_error(plot_ion_image(two, idx = 1, group = c("a", "b")), "one value per pixel")

  # wrappers pass `...` through
  expect_s3_class(
    plot_ion_image(two, idx = 1, group = FALSE, show_legend = FALSE),
    "ggplot"
  )
})
