test_that("default adducts contain required columns", {
  adducts <- default_adducts(NULL)

  expect_true(all(c("name", "mode", "mass", "sign") %in% names(adducts)))
  expect_true(all(adducts$mode %in% c("pos", "neg")))
  expect_true(is.numeric(adducts$mass))
  expect_true(is.numeric(adducts$sign))
})

test_that("default adducts reject invalid ion mode", {
  expect_error(default_adducts("bad_mode"))
})
