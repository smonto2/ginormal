#' Transformation of the PBDV function for use in the GIN distribution
#'
#' @param a degrees-of-freedom parameter.
#' @param z argument for the parabolyc cylinder function.
#' @param method string with the method used to compute the parabolic cylinder function.
#' method = "Fortran" uses the compiled Fortran version, while method = "R" uses an R translation.
#' @details
#' This function uses the PBDV Fortran subroutine provided in the
#' SPECFUN Fortran library by Shanjie Zhang and Jianming Jin in
#' Computation of Special Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
#'
#' @return Modified parabolyc cylinder function.
#'
#' @noRd
pbdam <- function(a, z, method = "Fortran") {
  am <- -a + 1
  return(pbdv_r(am, z, method))
}

