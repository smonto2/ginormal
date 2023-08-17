#' Apply mode shift for use in random variable generation.
#'
#' @param z quantile.
#' @param a degrees-of-freedom parameter.
#' @param m location parameter.
#' @param shift logical; should the argument be shifted by multiplying by
#' `z` minus the mode
#'
#' @return Shifted square root of the kernel.
#'
#' @noRd
mshift <- function(z, mode, alpha, mu, shift = TRUE){
  if (shift) {
    x <- z - mode
  }else{
    x <- 1
  }

  return(x * sqrt(dgin1(z, alpha, mu, FALSE, TRUE)))
}
