#Input:   x  --- Argument；va --- Order
#Output:  PV --- Vv(x)
#' Title
#'
#' @param va intermediate function's parameter
#' @param x intermediate function's parameter
#' @return vvla
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
