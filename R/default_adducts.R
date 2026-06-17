default_adducts <- function() {

  data.frame(
    name = c(
      # POSITIVE MODE
      "[M+H]+",
      "[M+Na]+",
      "[M+NH4]+",
      "[M+K]+",
      "[M+H-H2O]+",
      "[M+2Na-H]+",
      "[M+2K-H]+",

      # NEGATIVE MODE
      "[M-H]-",
      "[M+Cl]-",
      "[M-H-H2O]-"
    ),

    mode = c(
      rep("pos", 7),
      rep("neg", 3)
    ),

    # NET m/z shift (Da) — HARD-CODED
    mass = c(
      # pos
      +1.007276,    # [M+H]+
      +22.989218,  # [M+Na]+
      +18.033823,  # [M+NH4]+
      +38.963158,  # [M+K]+
      -17.003289,  # [M+H-H2O]+  (1.007276 - 18.010565)
      +44.971160,  # [M+2Na-H]+  (2*22.989218 - 1.007276)
      +76.919040,  # [M+2K-H]+   (2*38.963158 - 1.007276)

      # neg
      -1.007276,   # [M-H]-
      +34.968853,  # [M+Cl]-
      -19.017841   # [M-H-H2O]-  (1.007276 + 18.010565)
    ),

    # ALWAYS +1 (do not touch)
    sign = rep(+1, 10),

    stringsAsFactors = FALSE
  )
}


#To add a new adduct
#rbind(
#  default_adducts(),
#  data.frame(
#    name = "[M+FA-H]-",
#    mode = "neg",
#    mass = 43.9909,
#    sign = +1
#  )
#)
