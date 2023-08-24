#' Standardized generalized inverse normal density setting tau = 1
#' Internal function
#' @inheritParams dgin
#' @return Numeric scalar with density.
#' @noRd
dgin1 <- function(z, alpha, mu, log=TRUE, quasi=FALSE){
  am <- 0.5 * (alpha - 1)
  lkern <- (-alpha) * log(abs(z)) - 0.5 * (1 / z - mu) ^ 2
  if (quasi) {
    res <- lkern
  }else{
    hyp <- hypergeometric1F1(am, 0.5, 0.5*mu^2, log = TRUE)
    lcons <- -0.5 * (mu ^ 2) + am * log(2) + lgamma(am) + hyp
    res <- lkern - lcons
  }

  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

