# Generalized inverse normal (log) density
# Standardized density setting tau = 1
# Internal function
# This is an internal function that adds two numbers.
#' @param a  a degrees-of-freedom parameter
#' @param b A single numeric value.
#' @return dgin1
dgin1 <- function(z, a, m, log=TRUE, quasi=FALSE){
  if (z == 0 || a <= 1) {
    res = -Inf
  }else{
    am = 0.5 * (a - 1)
    lkern = (-a) * log(abs(z)) - 0.5 * (1 / z - m) ^ 2
    if (quasi==TRUE){
      res = lkern
    }else{
      hyp = hypergeometric1F1(am, 0.5, 0.5*m^2, log = TRUE)
      lcons = -0.5 * (m ^ 2)+ am * log(2)+lgamma(am)+hyp
      res = lkern - lcons
    }
  }

  if (log==TRUE){
    return(res)
  }else{
    return(exp(res))
  }
}







