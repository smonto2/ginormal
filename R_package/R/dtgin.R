#' Title
#'
#' @param z real number (Distribution is supported on real line, z can be any real number.
#' Currently, only scalar z is supported. Density is set to be 0 at z = 0.)
#' @param a  real number > 1 (Degrees of freedom parameter. Function returns 0 density if a <= 1.)
#' @param m  similar to a location parameter, it shifts the density of the distribution left and right
#' @param t similar to a scale parameter, it spreads the density of the distribution
#' @param sign sign = TRUE implies truncation to positive numbers and sign = FALSE to negative numbers
#' @param log should the logarithm of the density be returned? Defaults to TRUE.
#' @param quasi should the value of the kernel (or quasi-density) be returned? Defaults to FALSE
#'
#' @return dtgin
#' @export dtgin
dtgin <- function(z, a, m, t, sign=TRUE, log=TRUE, quasi=FALSE){
  if (a <= 1) {
    warning("alpha should be greater than 1")
  }
  if (t <= 0) {
    warning("tau should be greater than 0")
  }
  if (z < 0 && sign == TRUE) {
    warning("sign = True means truncation to positive reals; z should greater than 0")
  }
  if (z > 0 && sign == FALSE) {
    warning("sign = False means truncation to negative reals; z should less than 0")
  }

  if (sign==TRUE){
    t_z = (z <= 0)
    mult = -1
  }else{
    t_z = (z >= 0)
    mult = 1
  }

  if(t_z == TRUE || (a <= 1) || (t <= 0)){
    res = -Inf
  }else{
    mt = m / t
    lkern = dgin1(z * t, a, mt, TRUE, TRUE) + log(t)
    if (quasi==TRUE){
      res = lkern
    }else{
      lcons = -0.25 * mt ** 2 + log(gamma(a - 1)) + log(pbdam(a, mult * mt))
      res = lkern - lcons
    }
  }

  if (log==TRUE){
    return(res)
  }else{
    return(exp(res))
  }
}
