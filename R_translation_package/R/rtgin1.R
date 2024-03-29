#' Generating random numbers from the standardized generalized inverse normal distribution truncated to the positive or negative reals
#'
#' @inheritParams dtgin
#' @param size number of desired draws. Output is numpy vector of length equal to size)
#' @param algo string with desired algorithm to compute minimal bounding rectangle.
#' If "hormann", use the method from Hörmann and Leydold (2014). When "leydold", use the one from Leydold (2001).
#' Defaults to "hormann" and returns an error for any other values.
#' @param verbose logical; should the acceptance rate from the ratio-of-uniforms
#' method be provided along with additional information? Defaults to FALSE.
#' @details
#' Currently, only values of alpha > 2 are supported. For Bayesian posterior sampling,
#' alpha is always larger than 2 even for non-informative priors.
#' Generate from positive region (`z` > 0) hen `sign = TRUE`, and from
#' negative region (`z` < 0) when `sign = FALSE`. When `verbose = TRUE`,
#' a list is returned containing the actual draw in `value`, as well as average
#' acceptance rate `avg_arate` and total number of acceptance-rejection steps `ARiters`.
#'
#' @return If `verbose = FALSE` (default), a numeric vector of length `size`.
#' Otherwise, a list with components `value`, `avg_arate`, and `ARiters`
#' @export rtgin
#'
#' @noRd
rtgin1 <- function(a, m, sign, algo = 'hormann', verbose = FALSE){
  # If sign = True, draw from the positive region Z > 0
  # If sign = False, draw from the negative region Z < 0
  # Compute necessary values for ratio-of-uniforms method
  k <- sqrt(m^2 + 4 * a)
  if (sign) {
    mode <- (-m + k) / (2 * a)
    mult <- -1
  } else {
    mode <- (-m - k) / (2 * a)
    mult <- 1
  }

  if (algo == 'hormann') {
    algo_1 <- TRUE
  } else if (algo == 'leydold') {
    algo_1 <- FALSE
  } else {
    stop("algo must be either 'hormann' or 'leydold'")
  }

  # Draw from minimal bounding rectangle and accept draw
  vmax <- mshift(mode, mode, a, m, FALSE)
  if (algo_1) {
    # Hörmann and Leydold (2014) using Cardano method
    R <- polyroot(c(-mode, 1 + m * mode, -m + a * mode, 2 - a))
    roots <- Re(R)   #real part
    roots <- roots[(-mult*roots) > 0]
    uvals <- sapply(roots, mshift, mode = mode, a = a, m = m, shift = TRUE)
    uvals <- sort(uvals)
  } else {
    # Leydold (2001) using the proportionality constant
    lcons <- -0.25 * m ^ 2 + lgamma(a - 1) + log(pbdam(a, mult * m))
    vp <- exp(lcons) / vmax
    uvals <- c(-vp, vp)
  }
  # Acceptance-Rejection algorithm (ratio-of-uniforms)
  test <- FALSE
  counter <- 0
  max_iter <- 100
  while (!test) {
    u <- runif(1, uvals[1], uvals[2])
    v <- runif(1, 0, vmax)
    x <- (u / v) + mode
    test <- (2 * log(v) <= dgin1(x, a, m, TRUE, TRUE)) && (-mult * x > 0)
    counter <- counter + 1
    if (counter > max_iter) {
      x <- mode
      warning("No candidate draw found for these parameter values. Giving up an returning the mode of the distribution")
      break
    }
  }

  if (verbose) {
    return(list('value' = x, 'ARiters' = counter))
  } else {
    return(x)
  }
}
