
library(dplyr)
library(tidyr)

# This script builds a lookup table of theoretical isotope intensity ratios
# for the EIPS (Element Isotope Plausibility Score) model.
#
# For each element included in EIPS (N, O, S, Cl, Br), we precompute the expected
# isotopic intensity ratio (I_M+Δ / I_M0) as a function of the number of atoms (n),
# based exclusively on natural isotope abundances and probabilistic isotope models.
#
# Only the dominant isotopic channel closest to the monoisotopic peak is considered
# (M+1 for N, M+2 for O, S, Cl and Br), consistent with the isotope morphology
# detection step.
#
# The resulting table is independent of molecular mass and is used at runtime
# to compare theoretical and experimental isotope ratios without requiring
# access to any molecular database.

elements <- c("N","O","S","Cl","Br")

element_k <- c(
  N  = 1L,
  O  = 2L,
  S  = 2L,
  Cl = 2L,
  Br = 2L
)

element_delta <- c(
  N  = 0.997035,
  O  = 2.004245,
  S  = 1.99580,
  Cl = 1.99705,
  Br = 1.99795
)

element_model <- c(
  N  = "linear",
  O  = "linear",
  S  = "linear",
  Cl = "halogen",
  Br = "halogen"
)

# Map element -> abundance key in isotope_abundance
abund_key <- function(el) {
  paste0(
    el,
    if (el == "Cl") "37" else
      if (el == "Br") "81" else
        if (el == "S")  "34" else
          if (el == "O")  "18" else "15"
  )
}

theoretical_ratio <- function(n, p, model) {
  stopifnot(model %in% c("linear","halogen"))
  if (model == "linear")  return(n * p)
  # halogen channel: P(1)/P(0) = n p / (1-p)  (k=1, i.e., M+2)
  (n * p) / (1 - p)
}


# Extract formulas and compute n_max per element from DB

db <- ref_database %>%
  transmute(formula = MolecularFormula) %>%
  filter(!is.na(formula))

# Parse once per formula (faster than per-element parsing in loops)
# Store as named list (formula -> named vector of counts)
cf_list <- lapply(db$formula, parse_formula)

# Helper: get count for an element, always integer
get_n_el <- function(cf, el) {
  if (is.null(cf)) return(0L)

  # atomic named vector (numeric/integer with names)
  if (is.atomic(cf)) {
    nms <- names(cf)
    if (is.null(nms) || !length(nms) || !(el %in% nms)) return(0L)
    val <- cf[[el]]
    if (is.null(val) || is.na(val)) return(0L)
    return(as.integer(val))
  }

  # list-like objects (includes data.frame)
  if (is.list(cf)) {
    nms <- names(cf)
    if (is.null(nms) || !length(nms) || !(el %in% nms)) return(0L)
    val <- cf[[el]]  # safe now because we checked el %in% names(cf)
    if (is.null(val) || length(val) == 0 || all(is.na(val))) return(0L)
    return(as.integer(val[[1]]))
  }

  0L
}
# Compute n_max per element
n_max_df <- tibble(element = elements) %>%
  rowwise() %>%
  mutate(
    n_max = max(vapply(cf_list, get_n_el, integer(1), el = element), na.rm = TRUE)
  ) %>%
  ungroup()

# Optional safety cap (defendible) to avoid huge tables if DB contains outliers
# n_max_df <- n_max_df %>% mutate(n_max = pmin(n_max, 200L))


# Build eips_n_table (generic n -> R_theo)

eips_n_table <- n_max_df %>%
  mutate(
    model = element_model[element],
    k     = element_k[element],
    delta = element_delta[element],
    p     = vapply(element, function(el) isotope_abundance[[abund_key(el)]], numeric(1))
  ) %>%
  mutate(n = purrr::map(n_max, ~seq.int(1L, .x))) %>%
  tidyr::unnest(n) %>%
  mutate(
    R_theo = mapply(theoretical_ratio, n = n, p = p, model = model)
  ) %>%
  select(element, model, k, delta, p, n, R_theo) %>%
  arrange(element, n)


# Save to data/

usethis::use_data(eips_n_table, overwrite = TRUE)
