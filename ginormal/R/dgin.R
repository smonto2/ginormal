#' Density for the generalized inverse normal distribution
#'
#' @importFrom stats runif
#' @param z quantile.
#' @param alpha degrees-of-freedom parameter.
#' @param mu similar to location parameter, controls asymmetry of the distribution.
#' @param tau similar to scale parameter, controls spread of the distribution.
#' @param log logical; should the log of the density be returned? Defaults to TRUE.
#' @param quasi logical; should the quasi-density value be returned? Defaults to FALSE.
#' @details
#' Currently, only scalars are supported for the quantile and parameter values.
#' Density is supported on the entire real line, `z` and `mu` can take any value
#' in $(-\\infty, \\infty)$. Density is only defined for parameter values
#' `alpha`${} > 1$ or `tau`${} > 0$, so it is set to 0 outside of these values.
#' The quasi-density or kernel is the density without the normalization constant,
#' use `quasi = TRUE` for this behavior.
#'
#' @return Numeric scalar with density.
#' @export dgin
#'
#' @examples
#' # Computing (log) density
#' dgin(z = 1, alpha = 3, mu = 1, tau = 1, log = TRUE, quasi = FALSE)
#'
#' # Generalized inverse normal density with alpha = 5, mu = 0, tau = 1
#' z_vals <- seq(-5, 5, length.out = 200)
#' fz_unc <- sapply(z_vals, function(z) dgin(z, 5, 0, 1, FALSE))
#' plot(z_vals, fz_unc, type = "l", xlab = 'Values', ylab = 'Density')
dgin <- function(z, alpha, mu, tau, log=TRUE, quasi=FALSE){
  # Check parameter values
  if (alpha <= 1) {
    warning("alpha should be greater than 1")
    res <- -Inf
  } else if (tau <= 0) {
    warning("tau should be greater than 0")
    res <- -Inf
  } else if (z == 0) {
    res <- -Inf
  } else {
    # Compute density
    mt <- mu / tau
    if (quasi == TRUE){
      res <- dgin1(z * tau, alpha, mt, TRUE, TRUE) + alpha * log(tau)
    } else {
      res <- dgin1(z * tau, alpha, mt, TRUE, FALSE) + log(tau)
    }
  }

  # Return
  if (log == TRUE){
    return(res)
  } else {
    return(exp(res))
  }
}
