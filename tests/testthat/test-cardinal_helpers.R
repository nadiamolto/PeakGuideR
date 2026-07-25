test_that("cardinal_to_peakmatrix() works on a sparse-backend Cardinal object", {
  skip_if_not_installed("Cardinal")
  skip_if_not_installed("BiocParallel")

  suppressPackageStartupMessages(requireNamespace("Cardinal"))
  Cardinal::setCardinalBPPARAM(BiocParallel::SerialParam())

  set.seed(42)
  sim <- Cardinal::simulateImage(
    preset = 1, npeaks = 30, dim = c(5, 5), representation = "centroid"
  )
  sim <- Cardinal::normalize(sim, method = "tic")
  sim <- Cardinal::peakAlign(sim, tolerance = 5, units = "ppm")
  sim <- Cardinal::process(sim)

  # Real "processed"/peak-picked imzML data (as produced e.g. by an
  # upstream Xcalibur peak-picking step, converted through ProteoWizard)
  # ends up on this same sparse, on-disk-backed representation after
  # normalize()+peakAlign()+process() - not a synthetic edge case.
  intens_layer <- Cardinal::imageData(sim)[["intensity"]]
  expect_s4_class(intens_layer, "sparse_mat")

  pkm <- cardinal_to_peakmatrix(sim, dataset_name = "synthetic_sparse")

  expect_type(pkm, "list")
  expect_s3_class(pkm, "rMSIprocPeakMatrix")
  expect_true(is.matrix(pkm$intensity))
  expect_equal(length(pkm$mass), ncol(pkm$intensity))
  expect_equal(nrow(pkm$intensity), pkm$numPixels)
  expect_false(anyNA(pkm$intensity))
})
