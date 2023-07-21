#' Title Generator of generalized inverse normal variables
#'
#' @param size integer (Number of desired draws. Output is numpy vector of length equal to size.)
#' @param a a degrees-of-freedom parameter,
#' @param m similar to a location parameter, it shifts the density of the distribution left and right
#' @param t similar to a scale parameter, it spreads the density of the distribution
#' @param algo They take an additional argument algo, which can be either "hormann" or "leydold", and defaults to "hormann"
#' @return rgin
#' @export rgin
#'
#' @examples
#' rgin(1000, 5, 0, 1, algo=TRUE)
rgin <- function(size, a, m, t, algo='hormann'){
  if(a <= 2) {
    stop("alpha should be greater than 2")}
  if(t <= 0) {
    stop("tau should be greater than 0")}
  if(algo == 'hormann') {
    algo_1 <- TRUE
  } else if(algo == 'leydold') {
    algo_1 <- FALSE
  } else {
    stop("algo should be either 'hormann' or 'leydold'")
  }
  sign = (m >= 0)
  F0r = F0gin(a, m, t, sign)
  res = rep(0,size)
  side = (runif(size) <= F0r)
  sum_side = sum(side)
  res[side == 0] = rtgin(size - sum_side, a, m, t, !sign, algo)
  res[side == 1] = rtgin(sum_side, a, m, t, sign, algo)
  return(res)
}
