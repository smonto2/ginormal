#' Parabolic cylinder function for large arguments
#'
#' @param va order.
#' @param x argument.
#' @return Scalar with parabolic cylinder function of order va evaluated at x.
#' @details
#' This is an R translation of the DVLA Fortran subroutine provided in the
#' SPECFUN Fortran library by Shanjie Zhang and Jianming Jin in
#' Computation of Special Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
#'
#' @noRd
dvla <- function(va, x) {
  pi <- 3.141592653589793
  eps <- 1.0e-12
  a0 <- abs(x)^va * exp(-0.25 * x^2)
  r <- 1.0
  pd <- 1.0
  for (k in 1:16) {
    r <- -0.5 * r * (2.0*k - va - 1.0) * (2.0*k - va - 2.0) / (k * x^2)
    pd <- pd + r

    if (abs(r/pd) < eps) {
      break
    }
  }

  pd <- a0 * pd

  if (x < 0.0) {
    vl <- vvla(va, -x)
    pd <- pi * vl / gamma(-va) + cos(pi * va) * pd
  }

  return(pd)
}
