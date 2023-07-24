#Input:   x  --- Argument；va --- Order
#Output:  PD --- Dv(x)

#' Title
#'
#' @param va intermediate parameters
#' @param x intermediate parameters
#'
#' @return dvsa
dvsa <- function(va,x) {
  eps <- 1.0e-15
  sq2 <- sqrt(2.0)
  ep <- exp(-0.25 * x^2)
  va0 <- 0.5 * (1.0 - va)

  if (va == 0.0) {
    pd <- ep
  } else {
    if (x == 0.0) {
      if (va0 <= 0.0 && va0 == as.integer(va0)) {
        pd <- 0.0
      } else {
        ga0 <- gamma(va0)
        pd <- sqrt(pi) / (2.0^(-0.5 * va) * gamma(va0))
      }
    } else {
      g1 <- gamma(-va)
      a0 <- 2.0^(-0.5 * va - 1.0) * ep / g1
      pd <- gamma(-0.5 * va)
      r <- 1.0

      for (m in 1:250) {
        gm <- gamma(0.5 * (m - va))
        r <- -r * sq2 * x / m
        pd <- pd + gm * r
        if (abs(gm * r) < abs(pd) * eps) {
          break
        }
      }
      pd <- a0 * pd
    }
  }

  return(pd)
}
