library(dplyr)


#This code is intended to get the theoretical model for ech isotope considering the relative abundance
#which can be found in the isotope_constants.R file


element_model <- c(
  N  = "linear",
  O  = "linear",
  S  = "linear",
  Cl = "halogen",
  Br = "halogen"
)

# exact mass shifts used in iso_morphology_candidates (z=1)
element_delta <- c(
  N  = 0.997035,
  O  = 2.004245,
  S  = 1.99580,
  Cl = 1.99705,
  Br = 1.99795
)

element_k <- c(
  N  = 1L,
  O  = 2L,
  S  = 2L,
  Cl = 2L,
  Br = 2L
)

theoretical_ratio <- function(n, p, model) {
  stopifnot(all(model %in% c("linear", "halogen")))
  out <- numeric(length(n))

  idx_linear  <- model == "linear"
  idx_halogen <- model == "halogen"

  out[idx_linear]  <- n[idx_linear] * p
  out[idx_halogen] <- (n[idx_halogen] * p) / (1 - p)

  out
}

db <- ref_database %>%
  transmute(
    formula = MolecularFormula,
    mz_mono = MonoisotopicMass
  ) %>%
  filter(!is.na(formula), !is.na(mz_mono))

elements <- c("N", "O", "S", "Cl", "Br")

eips_table <- lapply(elements, function(el) {

  p <- isotope_abundance[paste0(
    el,
    if (el == "Cl") "37" else
      if (el == "Br") "81" else
        if (el == "S")  "34" else
          if (el == "O")  "18" else "15"
  )]

  model <- element_model[[el]]

  db %>%
    mutate(
      n_el = vapply(formula, function(f) {
        cf <- parse_formula(f)
        val <- if (el %in% names(cf)) cf[[el]] else 0
        val <- if (is.na(val)) 0 else val
        as.integer(val)
      }, integer(1)),            # <- integer
      element = el,
      model   = model,
      k       = element_k[[el]],      # <- explicit M+1/M+2
      delta   = element_delta[[el]],  # <- explicit Da shift
      R_theo  = theoretical_ratio(n_el, p, model)
    ) %>%
    filter(n_el > 0) %>%
    select(formula, mz_mono, element, model, k, delta, n_el, R_theo)
}) %>%
  bind_rows()

usethis::use_data(eips_table, overwrite = TRUE)
