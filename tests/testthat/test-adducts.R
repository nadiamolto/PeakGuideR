test_that("default adducts return expected ion modes", {
  pos <- default_adducts("pos")
  neg <- default_adducts("neg")
  all_adducts <- default_adducts(NULL)

  expect_true(is.data.frame(pos))
  expect_true(is.data.frame(neg))
  expect_true(is.data.frame(all_adducts))

  expect_true(all(pos$mode == "pos"))
  expect_true(all(neg$mode == "neg"))

  expect_true(all(c("name", "mode", "mass", "sign") %in% names(pos)))
  expect_true(all(c("name", "mode", "mass", "sign") %in% names(neg)))

  expect_gt(nrow(pos), 0)
  expect_gt(nrow(neg), 0)
})
