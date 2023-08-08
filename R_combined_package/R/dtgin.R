#' Density for the generalized inverse normal distribution truncated to the positive or negative reals
#'
#' @inheritParams dgin
#' @param sign logical. `TRUE` implies truncation to positive numbers (`z` > 0)
#' and `FALSE` to negative numbers (`z` < 0)
#' @param method two ways to calculate pbdv
#' @details
#' Currently, only scalars are supported for the quantile and parameter values.
#' Density is supported on the positive reals (`z` > 0) when `sign = TRUE` and to
#' negative reals (`z` < 0) when `sign = FALSE`. `mu` can take any value
#' in $(-\\infty, \\infty)$. Density is only defined for parameter values
#' `alpha`${} > 1$ or `tau`${} > 0$, so it is set to 0 outside of these values.
#' The quasi-density or kernel is the density without the normalization constant,
#' use `quasi = TRUE` for this behavior.
#'
#' @return Numeric scalar with density.
#' @export dtgin
dtgin <- function(z, alpha, mu, tau, sign=TRUE, log=TRUE, quasi=FALSE,method='Fortran'){
  # Check parameter values
  if (alpha <= 1) {
    warning("alpha should be greater than 1")
    res <- -Inf
  } else if (tau <= 0) {
    warning("tau should be greater than 0")
    res <- -Inf
  } else if ((z < 0) && (sign == TRUE)) {
    warning("sign = TRUE means truncation to positive reals; z should greater than 0")
    res <- -Inf
  } else if ((z > 0) && (sign == FALSE)) {
    warning("sign = FALSE means truncation to negative reals; z should less than 0")
    res <- -Inf
  } else if (z == 0) {
    res <- -Inf
  } else {
    # Compute truncated density
    mt = mu / tau
    lkern = dgin1(z * tau, alpha, mt, TRUE, TRUE) + log(tau)
    if (quasi == TRUE) {
      res <- lkern
    } else {
      mult <- 1 - 2*sign
      lcons <- -0.25 * mt ^ 2 + lgamma(alpha - 1) + log(pbdam(alpha, mult * mt,method))
      res <- lkern - lcons
    }
  }

  # Return
  if (log == TRUE){
    return(res)
  } else {
    return(exp(res))
  }
}
