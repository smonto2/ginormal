# Truncated generalized inverse normal log-densities
library(BAS)
#' Internal function
#'
#' This is an internal function that adds two numbers.
#'
#' @param a  a degrees-of-freedom parameter
#' @param z A single numeric value.
#' @return pbdam
pbdam <- function(a, z) {
  am = -a + 1
  return(pbdv_r(am,z))
}

