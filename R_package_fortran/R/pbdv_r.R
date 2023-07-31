#' title
#' @useDynLib GinFortran, .registration = TRUE
#' @param v ij
#' @param x ii
#'
#' @return pbdv_r(v,x)
pbdv_r <- function(v, x) {
  v <- as.double(v)
  x <- as.double(x)
  dv <- as.double(rep(0, 100))  # initialize an output array
  dp <- as.double(rep(0, 100))  # initialize an output array
  pdf <- as.double(0)  # initialize an output value
  pdd <- as.double(0)  # initialize an output value
  result <- .Fortran("pbdv", v = v, x = x, dv = dv, dp = dp, pdf = pdf, pdd = pdd, PACKAGE="GinFortran")

  return(result$pdf[length(result$pdf)])  # return the last value of pdf
}




