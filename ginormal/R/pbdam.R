#' Transformation of the PBDV function for use in the GIN distribution
#'
#' @param a degrees-of-freedom parameter.
#' @param z argument for the parabolyc cylinder function.
#' @param method For the method="R":This is an R translation of the PBDV
#' Fortran subroutine, provided in the SPECFUN Fortran library by Shanjie Zhang
#' and Jianming Jin in Computation of Special Functions, Wiley, 1996, ISBN:
#' -471-11963-6, LC:QA351.C45.
#' For the method="Fortran":This is the PBDV Fortran provided in the SPECFUN
#' Fortran library by Shanjie Zhang and Jianming Jin in
#' Computation of Special Functions, Wiley, 1996, ISBN: 0-471-11963-6, LC:
#' QA351.C45.
#' @return Modified parabolyc cylinder function.
#'
#' @noRd
pbdam <- function(a, z,method) {
  am = -a + 1
  return(pbdv_r(am,z,method))
}

