#' Generalized inverse normal (log) density
#' Standardized density setting tau = 1
#' Internal function
#' @inheritParams dgin
#' @return dgin1
dgin1 <- function(z, alpha, mu, log=TRUE, quasi=FALSE){
  if (z == 0 || alpha <= 1) {
    res = -Inf
  }else{
    am = 0.5 * (alpha - 1)
    lkern = (-alpha) * log(abs(z)) - 0.5 * (1 / z - mu) ^ 2
    if (quasi==TRUE){
      res = lkern
    }else{
      hyp = hypergeometric1F1(am, 0.5, 0.5*mu^2, log = TRUE)
      lcons = -0.5 * (mu ^ 2)+ am * log(2)+lgamma(am)+hyp
      res = lkern - lcons
    }
  }

  if (log==TRUE){
    return(res)
  }else{
    return(exp(res))
  }
}

