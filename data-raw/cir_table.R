
#This code is intended to build the CIR table:
#The main idea is to establish a theoretical M+0/M+1 ratio model for molecules with carbons (only C isotopes detected in the previous step will be evaluated with this method)

#CIR model bases:
# According to a binomial distribution:
# - Each C in each formula has an estimated probability to be 13C of 0.0107.
# - The probability that each  molecule has only one 13C (being the others 12C) is
# nC * P13C(1-P13C)^nC-1 --> nC*P13C when P13C is low (ref: Habler,K., Rexhaj, A., Adling-Ehrhardt, M., & Vogeser,M. (2024). Understanding isotopes, isomers, and isomars in mass spectrometry. Journal of Mass Spectrometry and Advances in the Clinical Lab, 33, 49-54.)
#This term is considered to be the theoretical I(M+1)/I(M+0) as is (according to the binomial distribution theory) aprox. P(1 13C)/P(0 13C)
#that is the same as nC * P13C(1-P13C)^nC-1 --> nC*P13C so I(M+0)/I(M+1) --> nC*P13C
# Cause: P(1)/P(0) = (nC*P13C)/ (1-P13C) --> nC*P13C

library(dplyr)

#Filter the fields of interest (formula and monoisotopic mass from my db)
db<- ref_database %>% transmute(formula= MolecularFormula,
                                mz_mono= MonoisotopicMass) %>%
  filter(!is.na(formula), !is.na(mz_mono))

#Count carbons in each formula
parse_formula <- function(formula) {
  m <- gregexpr("([A-Z][a-z]?)([0-9]*)", formula, perl = TRUE)
  parts <- regmatches(formula, m)[[1]]
  elems <- gsub("([A-Z][a-z]?)([0-9]*)", "\\1", parts, perl = TRUE)
  nums  <- gsub("([A-Z][a-z]?)([0-9]*)", "\\2", parts, perl = TRUE)
  nums[nums == ""] <- "1"
  counts <- as.numeric(nums)
  tapply(counts, elems, sum)
}

#Apply function getting an object with the number of C
nC <- vapply(db$formula, function(f) {
  cf <- parse_formula(f)
  if ("C" %in% names(cf)) cf[["C"]] else 0L
}, numeric(1))


# theoretical ratio mz_theo
db_C <- db %>%
  dplyr::mutate(
    nC = vapply(formula, function(f) {
      cf <- parse_formula(f)          # tu función de antes
      if ("C" %in% names(cf)) cf[["C"]] else 0L
    }, numeric(1)),
    R_theo = nC * 0.0107              # p_13C
  ) %>%
  filter(nC > 0)

# GAM model

library(mgcv)

fit_cir <- gam(
  R_theo ~ s(mz_mono, k = 20), #soft function
  #Different k were tested (30, 50, 80) and the lowest which allowed a stable curvature without artifacts was chosen.
  data   = db_C,
  family = gaussian(), #R theo is a continuous ratio, it is neither discrete nor proportional, so gaussian model is the most accurate.
  gamma=1.2 # To ensure a conservative curve (it penalizes if the curve is not soft)
)

gam.check(fit_cir)  # check if the parametres are fitted
#I can see that there is more variability for higher mz values
plot(fit_cir, shade = TRUE)  # plot


grid_mz <- seq(100, 1200, by = 1)  # it creates a vector of 1101 values.

cir_table <- data.frame(
  mz    = grid_mz,
  R_theo = predict(fit_cir,
                   newdata = data.frame(mz_mono = grid_mz),
                   type = "response")
)

#To get a teorethical ratio we have to interpolate values:
#approx(cir_table$mz, cir_table$R_theo, xout = 103.4567, rule = 2)$y
#approx: search the two nearest points to the selected value (in this case 103.4567) and calculates the relative position in the vector
#Then linealy interpolates the rations (theoretical vs estimated)
#rule=2 means that it will extrapolate for values not included in the model (ex: lower than 100 or higher than 1200)


usethis::use_data(cir_table, overwrite = TRUE)
