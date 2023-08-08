# Generator of truncated generalized inverse normal variables
#' Internal function
#'
#' This is an internal function that adds two numbers.
#'
#' @param z quantile
#' @param mode ?
#' @param m similar to location parameter, controls asymmetry of the distribution.
#' @param shift ?
#' @param a  a degrees-of-freedom parameter
#'
#' @return mshift
mshift <- function(z, mode, a, m, shift=TRUE){
  if (shift){
    x = z - mode
  }else{
    x = 1
  }

  return(x * sqrt(dgin1(z, a, m, FALSE, TRUE)))
}
