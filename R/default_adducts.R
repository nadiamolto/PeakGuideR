#' Default adduct definitions
#'
#' @description
#' Returns the default adduct definitions used by PeakGuideR for the selected
#' ion mode.
#'
#' @param ion_mode Ion mode. Either `"pos"` or `"neg"`. If `NULL`, adducts for
#'   both ion modes are returned.
#'
#' @return A data.frame containing adduct definitions, including adduct names,
#' ion mode and mass shifts relative to the neutral molecule.
#'
#' @examples
#' default_adducts("pos")
#' default_adducts("neg")
#' default_adducts(NULL)
#'
#' @export
default_adducts <- function(ion_mode = c("pos", "neg")) {

  adducts <- data.frame(
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

    # NET m/z shift (Da)
    mass = c(
      # pos
      +1.007276,   # [M+H]+
      +22.989218,  # [M+Na]+
      +18.033823,  # [M+NH4]+
      +38.963158,  # [M+K]+
      -17.003289,  # [M+H-H2O]+
      +44.971160,  # [M+2Na-H]+
      +76.919040,  # [M+2K-H]+

      # neg
      -1.007276,   # [M-H]-
      +34.968853,  # [M+Cl]-
      -19.017841   # [M-H-H2O]-
    ),

    # ALWAYS +1
    sign = rep(+1, 10),

    stringsAsFactors = FALSE
  )

  if (is.null(ion_mode)) {
    return(adducts)
  }

  ion_mode <- match.arg(ion_mode)

  adducts[adducts$mode == ion_mode, , drop = FALSE]
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
