#' Derivative of parabolic cylinder function in the notation of Whittaker
#'
#' @param v order.
#' @param x argument.
#' @return Scalar with derivative of parabolic cylinder function of order v evaluated at x.
#' @details
#' This is an R translation of the VVLA Fortran subroutine provided in the
#' SPECFUN Fortran library by Shanjie Zhang and Jianming Jin in
#' Computation of Special Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
#' Function can also produce derivatives of a given order, but this is not used.
#'
#' @noRd
vvla <- function(va, x) {
  eps <- 1.0e-12
  a0 <- abs(x)^(-va-1.0) * sqrt(2.0/pi) * exp(0.25 * x^2)
  r <- 1.0
  pv <- 1.0

  for (k in 1:18) {
    r <- 0.5 * r * (2.0*k + va - 1.0) * (2.0*k + va) / (k * x^2)
    pv <- pv + r

    if (abs(r/pv) < eps) {
      break
    }
  }
  pv <- a0 * pv

  if (x < 0.0) {
    pdl <- dvla(va, -x)
    dsl <- sin(pi * va)^2
    pv <- dsl * gamma(-va) / pi * pdl - cos(pi * va) * pv
  }

  return(pv)
}
