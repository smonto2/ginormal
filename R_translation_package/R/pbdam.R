#' Transformation of the PBDV function for use in the GIN distribution
#'
#' @param a degrees-of-freedom parameter.
#' @param z argument for the parabolyc cylinder function.
#' @return Modified parabolyc cylinder function.
#'
#' @noRd
library(BAS)
pbdam <- function(a, z) {
  am <- -a + 1
  return(pbdv(am,z))
}

