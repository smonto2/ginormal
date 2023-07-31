#' Generating random numbers from the generalized inverse normal distribution
#'
#' @inheritParams rtgin
#' @details
#' Currently, only values of alpha > 2 are supported. For Bayesian posterior sampling,
#' alpha is always larger than 2 even for non-informative priors. The algorithm requires
#' calculating the probability of truncation region (either z < 0 or z > 0). It is more stable
#' to compute a probability bounded away from 0. As mu controls asymmetry,
#' when mu > 0, P(truncation region) = P(z>0) >= 50%, and this probability is computed.
#' If mu < 0, P(z<0) >= 50% and this region's probability is used.
#'
#' @return Numeric vector of length `size`.
#' @export rgin
rgin <- function(size, alpha, mu, tau, algo = "hormann") {
    # Check parameter values (return error)
    if (alpha <= 2) {
        stop("alpha should be greater than 2")
    }
    if (tau <= 0) {
        stop("tau should be greater than 0")
    }
    if ((algo != "hormann") && (algo != "leydold")) {
        stop("algo should be either 'hormann' or 'leydold'")
    }

    # Generate using the algorithm for truncated variables
    sign <- (mu >= 0)
    F0r <- F0gin(alpha, mu, tau, sign)
    res <- rep(0, size)
    side <- (runif(size) <= F0r)
    sum_side <- sum(side)
    res[side == 0] <- rtgin(size - sum_side, alpha, mu, tau, !sign, algo)
    res[side == 1] <- rtgin(sum_side, alpha, mu, tau, sign, algo)
    return(res)
}
