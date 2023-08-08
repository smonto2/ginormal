#' Generating random numbers from the generalized inverse normal distribution truncated to the positive or negative reals
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
rtgin <- function(size, alpha, mu, tau, sign, algo='hormann', verbose=FALSE) {
  # Check parameter values (return error)
  if (alpha <= 2) {
    stop("alpha should be greater than 2")
  }
  if (tau <= 0) {
    stop("tau should be greater than 0")
  }
  if ((algo != 'hormann') && (algo != 'leydold')) {
    stop("algo should be either 'hormann' or 'leydold'")
  }

  if (size == 0) {
    # When no size, return a null
    return(NULL)
  } else {
    # Generate using standardized kernel
    res <- rep(0, size)
    if (verbose) {
      ARiters <- rep(0, size)
    }
    mt <- mu / tau
    for (i in 1:size) {
      temp <- rtgin1(alpha, mt, sign, algo, verbose)
      if (verbose) {
        res[i] <- temp$value / tau
        ARiters[i] <- temp$ARiters
      } else {
        res[i] <- temp / tau
      }
    }
    if (verbose) {
      return(list(value = res, avg_arate = mean(1/ARiters), ARiters = ARiters))
    } else {
      return(res)
    }
  }
}

