#' Title
#'
#' @param size integer
#' (Number of desired draws. Output is numpy vector of length equal to size)
#' @param a real number > 2
#' (Degrees of freedom parameter. Currently, only values of a > 2 are supported.m : real number)
#' @param m real number (Similar to location parameter, controls asymmetry of the distribution.)
#' @param t similar to a scale parameter, it spreads the density of the distribution
#' @param sign where sign = TRUE implies truncation to positive numbers and sign = FALSE to negative numbers
#' @param algo boolean, optional
#' algo = True, use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
#' algo = False, use the one from Leydold (2001)
#' Defaults to True.
#' @param verbose boolean, optional
#' Should the acceptance rate from the ratio-of-uniforms method be provided
#' along with additional information? Defaults to False.
#' @return rtgin
#' @export rtgin
rtgin <- function(size, a, m, t, sign, algo, verbose=FALSE) {
  if (a <= 2) {
    stop("alpha should be greater than 2")
  }
  if (t <= 0) {
    stop("tau should be greater than 0")
  }
  if (algo == 'hormann') {
    algo_1 <- TRUE
  } else if (algo == 'leydold') {
    algo_1 <- FALSE
  } else {
    stop("algo should be either 'hormann' or 'leydold'")
  }

  res <- rep(0, size)
  if (verbose) {
    ARiters <- rep(0, size)
  }
  mt <- m / t
  for (i in 1:size) {
    temp <- rtgin1(a, mt, sign, algo, verbose)
    if (verbose) {
      res[i] <- temp$value / t
      ARiters[i] <- temp$ARiters
    } else {
      res[i] <- temp / t
    }
  }
  if (verbose) {
    return(list(value=res, avg_arate=mean(1/ARiters), ARiters=ARiters))
  } else {
    return(res)
  }
}
