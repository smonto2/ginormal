#' Title General density with tau from (0, infinity)
#' @importFrom stats runif
#' @param z real number (Distribution is supported on real line, z can be any real number.Currently, only scalar z is supported. Density is set to be 0 at z = 0.)
#' @param a a degrees-of-freedom parameter
#' @param m similar to a location parameter, it shifts the density of the distribution left and right
#' @param t real number > 0 (Similar to scale parameter, controls spread of the distribution.)
#' @param log boolean, optional(Should the log of the density be returned? Default is True.)
#' @param quasi boolean, optional(The quasi-density or kernel is the density without the normalization constant.
#' Should the quasi-density value be returned? Default is False.)
#' @return dgin
#' @export dgin

dgin <- function(z, a, m, t, log=TRUE, quasi=FALSE){
  if (a <= 1) {
    warning("alpha should be greater than 1")
  }
  if (t <= 0) {
    warning("tau should be greater than 0")
  }

  if (t<=0){
    res <- -Inf
  }else{
    mt <- m / t
    if (quasi==TRUE){
      res = dgin1(z * t, a, mt, TRUE, TRUE) + a * log(t)
    }else{
      res = dgin1(z * t, a, mt, TRUE, FALSE) + log(t)
    }
  }

  if (log==TRUE){
    return(res)
  }else{
    return(exp(res))
  }
}
