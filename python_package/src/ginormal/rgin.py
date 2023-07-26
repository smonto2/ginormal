import numpy as np
from scipy import stats
from ginormal.utils.F0gin import F0gin
from ginormal.rtgin import rtgin

def rgin(size, alpha, mu, tau, algo='hormann'):
    """
    Generating random numbers from the generalized inverse normal distribution
    supported on the real numbers. Currently, only values of alpha > 2 are supported.
    For Bayesian posterior sampling, alpha is always larger than 2 even for non-informative priors.
    As mu controls asymmetry, when mu > 0, P(truncation region) = P(z>0) >= 50%,
    and this probability is computed. If mu < 0, P(z<0) >= 50% and this region's probability is used. 

    Parameters
    ----------
    size : integer
        Number of desired draws. Output is numpy vector of length equal to size.
    alpha : real number > 2
        Degrees of freedom parameter. Currently, only values of alpha > 2 are supported.
    mu : real number
        Similar to location parameter, controls asymmetry of the distribution.
    tau : real number > 0
        Similar to scale parameter, controls spread of the distribution.
    algo : string, optional
        algo = 'hormann', use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
        algo = 'leydold', use the one from Leydold (2001)
        Defaults to 'hormann'.
    """
    if alpha <= 2:
        raise ValueError("alpha should be greater than 2")
    if tau <= 0:
        raise ValueError("tau should be greater than 0")
    if (algo != 'hormann') & (algo != 'leydold'):
        raise ValueError("algo should be either 'hormann' or 'leydold'")
    sign = (mu >= 0)
    F0r = F0gin(alpha, mu, tau, sign)
    res = np.zeros(size)
    side = (stats.uniform.rvs(size=size) <= F0r)
    sum_side = sum(side)
    res[side == 0] = rtgin(size - sum_side, alpha, mu, tau, not sign, algo)
    res[side == 1] = rtgin(sum_side, alpha, mu, tau, sign, algo)
    return res
