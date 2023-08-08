# CDF of generalized inverse normal at 0
#' Internal function
#' This is an internal function that adds two numbers.
#' @importFrom BAS hypergeometric1F1
#' @param a  a degrees-of-freedom parameter
#' @param m A single numeric value.
#' @param t parameter
#' @param method two ways to calculate pbdv
#' @param sign default TRUE
#' @return F0gin
F0gin <- function(a, m, t, sign=TRUE,method="Fortran"){
  mult = 1 - 2 * sign
  mt = m / t
  parhyp = log(pbdam(a, mult * mt,method)) - hypergeometric1F1(0.5 * (a - 1), 0.5, 0.5*mt^2, log = TRUE)
  lcons = 0.25 * (mt ^ 2) + 0.5 * (a - 3) * log(2) - 0.5 * log(pi) + lgamma(0.5 * a) + parhyp
  return(exp(lcons))
}
