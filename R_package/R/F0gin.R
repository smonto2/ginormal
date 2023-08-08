#' Probability of truncation region for truncated GIN distribution
#'
#' @importFrom BAS hypergeometric1F1
#' @param a degrees-of-freedom parameter.
#' @param m location parameter.
#' @param t shape parameter.
#' @param sign logical. `TRUE` implies truncation to positive numbers (`z` > 0)
#' and `FALSE` to negative numbers (`z` < 0)
#' @return F0gin
#'
#' @noRd
F0gin <- function(a, m, t, sign=TRUE){
  mult = 1 - 2 * sign
  mt = m / t
  parhyp = log(pbdam(a, mult * mt)) - hypergeometric1F1(0.5 * (a - 1), 0.5, 0.5*mt^2, log = TRUE)
  lcons = 0.25 * (mt ^ 2) + 0.5 * (a - 3) * log(2) - 0.5 * log(pi) + lgamma(0.5 * a) + parhyp
  return(exp(lcons))
}
