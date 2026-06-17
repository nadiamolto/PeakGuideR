#' Parse molecular formula into element counts
#'
#' @param formula Character molecular formula (e.g. "C6H12O6")
#'
#' @return Named numeric vector with element counts
#' @keywords internal
parse_formula <- function(formula) {
  m <- gregexpr("([A-Z][a-z]?)([0-9]*)", formula, perl = TRUE)
  parts <- regmatches(formula, m)[[1]]
  elems <- sub("([A-Z][a-z]?)([0-9]*)", "\\1", parts, perl = TRUE)
  nums  <- sub("([A-Z][a-z]?)([0-9]*)", "\\2", parts, perl = TRUE)
  nums[nums == ""] <- "1"
  tapply(as.numeric(nums), elems, sum)
}
