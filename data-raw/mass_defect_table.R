#This code is intended to build the mass defect reference table used to flag inferred neutral masses
#(see mass_defect_ok() in R/mass_defect.R).
#The main idea is to describe, for each integer mass, which mass defects are observed among known metabolites.

#Mass defect definition (classical): the decimal part of the neutral mass, mass - floor(mass), a value in [0, 1).
#The integer part is the integer mass the neutral mass belongs to (e.g. 130.0 and 130.9 both belong to 130).

#Reference range: for each integer mass from 55 to 1300, the minimum and the maximum mass defect observed among the HMDB
#and ChEBI entries of the compound mass database. Integer masses without any entry get NA limits (not evaluable).

library(dplyr)

devtools::load_all() #to access the internal reference tables already stored in R/sysdata.rda

#ref_database is the complete compound mass database (HMDB, ChEBI and NORMAN entries), as in the other data-raw scripts
db <- ref_database %>%
  filter(Source %in% c("HMDB", "CHEBI"), !is.na(MonoisotopicMass)) %>%
  transmute(mass = MonoisotopicMass)

grid_min <- 55
grid_max <- 1300

db <- db %>%
  mutate(integer_mass = floor(mass),
         defect = mass - integer_mass) %>%
  filter(integer_mass >= grid_min, integer_mass <= grid_max)

mass_defect_table <- data.frame(mass = grid_min:grid_max) %>%
  left_join(
    db %>%
      group_by(mass = integer_mass) %>%
      summarise(lower = min(defect), upper = max(defect), n = n(), .groups = "drop"),
    by = "mass"
  ) %>%
  mutate(n = ifelse(is.na(n), 0L, n))

#The internal reference tables share R/sysdata.rda, so all of them are saved together
cir_table <- get("cir_table", envir = asNamespace("PeakGuideR"))
eips_table <- get("eips_table", envir = asNamespace("PeakGuideR"))
eips_n_table <- get("eips_n_table", envir = asNamespace("PeakGuideR"))

usethis::use_data(cir_table, eips_table, eips_n_table, mass_defect_table,
                  internal = TRUE, overwrite = TRUE)
