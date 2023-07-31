path <- system.file(package="Gin")
#path <- getwd()
compile_fortran <- function() {
  fortran_file_path <- paste0(path, "./src/fortfuncsf.f90")
  shared_object_path <- paste0(path, "./src/fortfuncs.so")
  system(paste("gfortran -fpic -shared", fortran_file_path, "-o", shared_object_path))
}
# Call the function to compile the Fortran code
compile_fortran()
#setwd(path)
dyn.load("./src/fortfuncs.so")
#' Title
#'
#' @param v ij
#' @param x ii
#'
#' @return pbdv(v,x)
#' @export
pbdv <- function(v, x) {
  v <- as.double(v)
  x <- as.double(x)
  dv <- as.double(rep(0, 100))  # initialize an output array
  dp <- as.double(rep(0, 100))  # initialize an output array
  pdf <- as.double(0)  # initialize an output value
  pdd <- as.double(0)  # initialize an output value
  result <- .Fortran("pbdv", v = v, x = x, dv = dv, dp = dp, pdf = pdf, pdd = pdd, PACKAGE="GinFortran")

  return(result$pdf[length(result$pdf)])  # return the last value of pdf
}




