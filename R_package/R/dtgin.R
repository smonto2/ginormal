#' Density for the generalized inverse normal distribution truncated to the positive or negative reals
#'
#' @inheritParams dgin
#' @param sign logical. `TRUE` implies truncation to positive numbers (`z` > 0)
#' and `FALSE` to negative numbers (`z` < 0).
#' @param method string with the method used to compute the parabolic cylinder function
#' in the normalization constant. `method = "Fortran"` uses a compiled Fortran version,
#' which is the default. `method = "R"` uses an R translation of this function.
#' @details
#' Currently, only scalars are supported for the quantile and parameter values.
#' Density is supported on the positive reals (`z` > 0) when `sign = TRUE` and to
#' negative reals (`z` < 0) when `sign = FALSE`. `mu` can take any value
#' in \eqn{(-\infty, \infty)}. Density is only defined for parameter values
#' `alpha` > 1 or `tau` > 0, so it is set to 0 outside of these values.
#' The quasi-density or kernel is the density without the normalization constant,
#' use `quasi = TRUE` for this behavior.
#'
#' @return Numeric scalar with density.
#' @export dtgin
#'
#' @examples
#' # Computing (log) truncated densities
#' dtgin(z = 1, alpha = 3, mu = 1, tau = 1, sign = TRUE, log = TRUE, quasi = FALSE)
#' dtgin(z = -1, alpha = 3, mu = -1, tau = 1, sign = FALSE, log = TRUE, quasi = FALSE)
#'
#' # Generalized inverse normal density with alpha = 5, mu = 0, tau = 1
#' n_values <- 200
#' z_vals <- seq(-5, 5, length.out = n_values)
#'
#' # Truncated to positive reals (z > 0)
#' fz_p <- sapply(z_vals[z_vals > 0], function(z) dtgin(z, 5, 0, 1, TRUE, FALSE))
#' fz_p <- c(rep(0, n_values - sum(z_vals > 0)), fz_p)
#' plot(z_vals, fz_p, type = "l", xlab = 'Values', ylab = 'Density')
#'
#' # Truncated to positive reals (z < 0)
#' fz_n <- sapply(z_vals[z_vals < 0], function(z) dtgin(z, 5, 0, 1, FALSE, FALSE))
#' fz_n <- c(fz_n, rep(0, n_values - sum(z_vals < 0)))
#' plot(z_vals, fz_n, type = "l", xlab = 'Values', ylab = 'Density')
#'
#' # Both truncated densities together
#' plot(z_vals, fz_p, type = "l", xlab = 'Values', ylab = 'Density')
#' lines(z_vals, fz_n, col = 'blue', lty = 2)
#' legend('topright', legend = c('z > 0', 'z < 0'),
#'        col = c('black', 'blue'), lty = 1:2)
dtgin <- function(z, alpha, mu, tau, sign = TRUE, log = TRUE,
                  quasi = FALSE, method = 'Fortran'){
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
    mt <- mu / tau
    lkern = dgin1(z * tau, alpha, mt, TRUE, TRUE) + log(tau)
    if (quasi) {
      res <- lkern
    } else {
      mult <- 1 - 2*sign
      lcons <- -0.25 * mt ^ 2 + lgamma(alpha - 1) + log(pbdam(alpha, mult * mt, method))
      res <- lkern - lcons
    }
  }

  # Return
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}
