test_that("mass defect reference table has the expected structure", {
  tab <- get("mass_defect_table", envir = asNamespace("PeakGuideR"))

  expect_equal(names(tab), c("mass", "lower", "upper", "n"))
  expect_equal(tab$mass, 55:1300)

  with_ref <- tab$n > 0
  expect_true(any(with_ref))
  expect_true(all(tab$lower[with_ref] >= 0 & tab$upper[with_ref] < 1))
  expect_true(all(tab$lower[with_ref] <= tab$upper[with_ref]))

  # integer masses without reference compounds have no limits
  expect_true(all(is.na(tab$lower[!with_ref])))
  expect_true(all(is.na(tab$upper[!with_ref])))
  expect_true(all(!is.na(tab$lower[with_ref])))
})

test_that("mass_defect_value is the decimal part of the mass", {
  expect_equal(mass_defect_value(180.0634), 0.0634, tolerance = 1e-8)
  expect_equal(mass_defect_value(c(130, 130.9)), c(0, 0.9), tolerance = 1e-8)

  # always between 0 and 1, including masses just below an integer
  expect_equal(mass_defect_value(c(97.9769, 134.0215)), c(0.9769, 0.0215), tolerance = 1e-8)
  expect_equal(mass_defect_value(900.8146), 0.8146, tolerance = 1e-8)
  expect_true(all(mass_defect_value(c(55.5, 301.99, 1299.01)) >= 0 &
                    mass_defect_value(c(55.5, 301.99, 1299.01)) < 1))
})

test_that("mass_defect_limits uses the row of the integer mass", {
  tab <- get("mass_defect_table", envir = asNamespace("PeakGuideR"))

  # 180.06 and 180.99 belong to 180; 181.0 belongs to 181
  lim <- mass_defect_limits(c(180.06, 180.99, 181.0))
  expect_equal(lim$lower, tab$lower[match(c(180, 180, 181), tab$mass)])
  expect_equal(lim$upper, tab$upper[match(c(180, 180, 181), tab$mass)])

  # outside the table or missing
  out <- mass_defect_limits(c(54.9, 1301, NA))
  expect_true(all(is.na(out$lower)))
  expect_true(all(is.na(out$upper)))

  # an integer mass without reference compounds is not evaluable
  empty <- tab$mass[tab$n == 0]
  if (length(empty)) {
    no_ref <- mass_defect_limits(empty[1] + 0.5)
    expect_true(is.na(no_ref$lower))
    expect_true(is.na(no_ref$upper))
  }
})

test_that("mass_defect_ok flags accepted, not accepted and not evaluable masses", {
  tab <- get("mass_defect_table", envir = asNamespace("PeakGuideR"))

  # typical metabolites
  expect_true(all(mass_defect_ok(c(180.0634, 307.0838, 134.0215, 268.0808))))

  # the limits themselves are accepted; a defect above the upper limit is not
  row <- tab[which(tab$n > 0 & tab$upper < 0.9)[1], ]
  expect_true(mass_defect_ok(row$mass + row$lower))
  expect_true(mass_defect_ok(row$mass + row$upper))
  expect_false(mass_defect_ok(row$mass + row$upper + 0.05))

  # outside the table range or missing
  expect_true(all(is.na(mass_defect_ok(c(45, 1500, NA)))))
})

test_that("add_mass_defect_flag appends one column and leaves the rest untouched", {
  tab <- get("mass_defect_table", envir = asNamespace("PeakGuideR"))
  row <- tab[which(tab$n > 0 & tab$upper < 0.9)[1], ]
  outside <- row$mass + row$upper + 0.05

  x <- data.frame(id = 1:3, neutral_mass_consensus = c(180.0634, outside, NA))
  y <- add_mass_defect_flag(x)

  expect_equal(names(y), c("id", "neutral_mass_consensus", "mass_defect_ok"))
  expect_equal(y$mass_defect_ok, c(TRUE, FALSE, NA))
  expect_equal(y[, c("id", "neutral_mass_consensus")], x)

  # zero-row table
  z <- add_mass_defect_flag(x[0, ])
  expect_equal(nrow(z), 0L)
  expect_true("mass_defect_ok" %in% names(z))

  # tables without neutral_mass_consensus are returned unchanged
  expect_identical(add_mass_defect_flag(data.frame(a = 1)), data.frame(a = 1))
})

test_that("workflow adds mass_defect_ok without changing the existing outputs", {
  data("example_pkm", package = "PeakGuideR")

  res_on <- run_peakguider_workflow(
    pkm = example_pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE
  )
  res_off <- run_peakguider_workflow(
    pkm = example_pkm, ion_mode = "pos", matrix = "HCCA", quiet = TRUE,
    mass_defect = FALSE
  )

  for (tbl in c("neutral_mass_candidates", "candidate_annotations")) {
    on <- res_on[[tbl]]
    off <- res_off[[tbl]]

    expect_false("mass_defect_ok" %in% names(off))
    expect_equal(tail(names(on), 1), "mass_defect_ok")
    expect_type(on$mass_defect_ok, "logical")
    expect_equal(on[, names(off), drop = FALSE], off)
  }

  expect_true(res_on$parameters$mass_defect)
  expect_false(res_off$parameters$mass_defect)
  expect_false("mass_defect_band" %in% names(res_on$parameters))
})
