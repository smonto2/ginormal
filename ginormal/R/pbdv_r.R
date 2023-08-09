#' Parabolic cylinder function in the notation of Whittaker
#'
#' @useDynLib ginormal, .registration = TRUE
#' @param v order.
#' @param x argument.
#' @return Scalar with parabolic cylinder function of order v evaluated at x.
#' @param method For the method="R":This is an R translation of the PBDV
#' Fortran subroutine, provided in the SPECFUN Fortran library by Shanjie Zhang
#' and Jianming Jin in Computation of Special Functions, Wiley, 1996, ISBN:
#' -471-11963-6, LC:QA351.C45.
#' For the method="Fortran":This is the PBDV Fortran provided in the SPECFUN
#' Fortran library by Shanjie Zhang and Jianming Jin in
#' Computation of Special Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC:
#' QA351.C45.
#' Function can also produce derivatives of a given order, but this is not used
#' .
#' @noRd
pbdv_r <- function(v, x, method = "Fortran"){
  if (method == "Fortran") {
    v <- as.double(v)
    x <- as.double(x)
    dv <- as.double(rep(0, 100))  # initialize an output array
    dp <- as.double(rep(0, 100))  # initialize an output array
    pdf <- as.double(0)  # initialize an output value
    pdd <- as.double(0)  # initialize an output value
    result <- .Fortran("pbdv", v = v, x = x, dv = dv, dp = dp, pdf = pdf, pdd = pdd, PACKAGE="ginormal")
    return(result$pdf[length(result$pdf)])  # return the last value of pdf
  } else if (method == "R") {
    xa <- abs(x)
    dv <- as.numeric(list())
    dp <- as.numeric(list())
    vh <- v
    v <- v + sign(v)
    nv <- as.integer(v)
    v0 <- v - nv
    na<- abs(nv)
    ep <- exp(-0.25 * x^2)
    ja <- 0
    if ( na>= 1) {
      ja <- 1
    }

    if (v >= 0.0) {
      if (v0 == 0.0) {
        pd0 <- ep
        pd1 <- x * ep
      } else {
        for (l in 0:ja) {
          v1 <- v0 + l
          if (xa <= 5.8) {
            pd1 <- dvsa(v1, x)
          } else {
            pd1 <- dvla(v1, x)
          }
          if (l == 0) {
            pd0 <- pd1
          }
        }
      }

      dv[0] <- pd0
      dv[1] <- pd1

      for (k in 2:na) {
        pdf <- x * pd1 - (k + v0 - 1.0) * pd0
        dv[k] <- pdf
        pd0 <- pd1
        pd1 <- pdf
      }
    } else {
      if (x <= 0.0) {
        if (xa <= 5.8) {
          pd0 <- dvsa(v0, x)
          v1 <- v0 - 1.0
          pd1 <- dvsa(v1, x)
        } else {
          pd0 <- dvla(v0, x)
          v1 <- v0 - 1.0
          pd1 <- dvla(v1, x)
        }
        dv[0] <- pd0
        dv[1] <- pd1

        for (k in 2:na) {
          pd <- (-x * pd1 + pd0) / (k - 1.0 - v0)
          dv[k] <- pd
          pd0 <- pd1
          pd1 <- pd
        }
      } else if (x <= 2.0) {
        v2 <- nv + v0
        if (nv == 0) {
          v2 <- v2 - 1.0
        }
        nk <- as.integer(-v2)
        f1 = dvsa(v2, x)
        v1 <- v2 + 1.0
        f0 <- dvsa(v1, x)
        dv[nk] <- f1
        dv[nk - 1] <- f0

        for (k in nk-2:0) {
          f <- x * f0 + (k - v0 + 1.0) * f1
          dv[k] <- f
          f1 <- f0
          f0 <- f
        }
      } else {
        if (xa <= 5.8) {
          pd0 <- dvsa(v0, x)
        } else {
          pd0 <- dvla(v0, x)
        }
        dv[0] <- pd0
        m <- 100 + na
        f1 <- 0.0
        f0 <- 1.0e-30
        f <- 0.0

        for (k in m:0) {
          f <- x * f0 + (k - v0 + 1.0) * f1
          if (k <=na ) {
            dv[k] <- f
          }
          f1 <- f0
          f0 <- f
        }

        s0 <- pd0 / f

        for (k in 0:na) {
          dv[k] <- s0 * dv[k]
        }
      }
    }
    for (k in 0:(na-1)) {
      v1 <- abs(v0) + k
      if (v >= 0.0) {
        dp[k] <- 0.5 * x * dv[k] - dv[k+1]
      } else {
        dp[k] <- -0.5 * x * dv[k] - v1 * dv[k+1]
      }
    }

    pdf <- dv[-1]
    pdd <- dp[-1]
    v <- vh
    result <- list(pdf)
    return(result[[1]][length(result[[1]]) - 1])
  } else {
    stop("Unknown method for computing parabolic cylinder function.")
  }
}
