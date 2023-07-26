import numpy as np
import warnings
from scipy import special, stats
from ginormal.utils.dgin1 import dgin1
from ginormal.utils.rtgin1 import rtgin1
from ginormal.utils.F0gin import F0gin
from ginormal.utils.pbdam import pbdam

def dgin(z, alpha, mu, tau, log=True, quasi=False):
    """
    Density for the generalized inverse normal distribution.

    Parameters
    ----------
    z : real number
        Distribution is supported on real line, z can be any real number.
        Currently, only scalar z is supported. Density is set to be 0 at z = 0.
    alpha : real number > 1
        Degrees of freedom parameter. Function returns 0 density if alpha <= 1.
    mu : real number
        Similar to location parameter, controls asymmetry of the distribution.
    tau : real number > 0
        Similar to scale parameter, controls spread of the distribution.
    log : boolean, optional
        Should the log of the density be returned? Default is True.
    quasi : boolean, optional
        The quasi-density or kernel is the density without the normalization constant.
        Should the quasi-density value be returned? Default is False.
    """
    # Check parameter values
    if alpha <= 1:
        warnings.warn("alpha should be greater than 1")
        res = -np.inf
    elif tau <= 0:
        warnings.warn("tau should be greater than 0")
        res = -np.inf
    elif z == 0:
        res = -np.inf
    else:
        # Compute density
        mt = mu / tau
        if quasi:
            res = dgin1(z * tau, alpha, mt, True, True) + alpha * np.log(tau)
        else:
            res = dgin1(z * tau, alpha, mt, True, False) + np.log(tau)
    
    # Return
    if log:
        return res
    else:
        return np.exp(res)

def dtgin(z, alpha, mu, tau, sign=True, log=True, quasi=False):    
    """
    Density for the generalized inverse normal distribution truncated to
    the positive or negative reals.

    Parameters
    ----------
    z : real number
        Distribution is supported on real line, z can be any real number.
        Currently, only scalar z is supported. Density is set to be 0 at z = 0.
    alpha : real number > 1
        Degrees of freedom parameter. Function returns 0 density if alpha <= 1.
    mu : real number
        Similar to location parameter, controls asymmetry of the distribution.
    tau : real number > 0
        Similar to scale parameter, controls spread of the distribution.
    sign : boolean
        sign = True means to truncate to positive reals (z > 0)
        sign = False truncates to negative reals (z < 0)
    log : boolean, optional
        Should the log of the density be returned? Default is True.
    quasi : boolean, optional
        The quasi-density or kernel is the density without the normalization constant.
        Should the quasi-density value be returned? Default is False.
    """
    # Check parameter values
    if alpha <= 1:
        warnings.warn("alpha should be greater than 1")
        res = -np.inf
    elif tau <= 0:
        warnings.warn("tau should be greater than 0")
        res = -np.inf
    elif (z < 0) and (sign == True):
        warnings.warn("sign = True means truncation to positive reals; z should be greater than 0")
        res = -np.inf
    elif (z > 0) and (sign == False):
        warnings.warn("sign = False means truncation to negative reals; z should be less than 0")
        res = -np.inf
    elif z == 0:
        res = -np.inf
    else:
        # Compute truncated density
        mt = mu / tau
        lkern = dgin1(z * tau, alpha, mt, True, True) + np.log(tau)
        if quasi:
            res = lkern
        else:
            mult = 1 - 2*sign
            lcons = -0.25 * mt ** 2 + special.loggamma(alpha - 1) + np.log(pbdam(alpha, mult * mt))
            res = lkern - lcons
    
    # Return
    if log:
        return res
    else:
        return np.exp(res)

def rtgin(size, alpha, mu, tau, sign=True, algo='hormann', verbose=False):
    """
    Generating random numbers from the generalized inverse normal distribution
    truncated to the positive or negative reals. Currently, only values of alpha > 2 are supported.
    For Bayesian posterior sampling, alpha is always larger than 2 even for non-informative priors.

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
    sign : boolean
        sign = True means to draw from positive reals, imposing z > 0
        sign = False draws from negative reals, z < 0
    algo : string, optional
        algo = 'hormann', use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
        algo = 'leydold', use the one from Leydold (2001)
        Defaults to 'hormann'.
    verbose : boolean, optional
        Should the acceptance rate from the ratio-of-uniforms method be provided
        along with additional information? Defaults to False. 
    """
    # Check parameter values (return error)
    if alpha <= 2:
        raise ValueError("alpha should be greater than 2")
    if tau <= 0:
        raise ValueError("tau should be greater than 0")
    if (algo != 'hormann') & (algo != 'leydold'):
        raise ValueError("algo should be either 'hormann' or 'leydold'")
    
    # Generate using standardized kernel
    res = np.zeros(size)
    if verbose: ARiters = np.copy(res)
    mt = mu / tau
    for i in range(size):
        temp = rtgin1(alpha, mt, sign, algo, verbose)
        if verbose:
            res[i] = temp['value'] / tau
            ARiters[i] = temp['ARiters']
        else:
            res[i] = temp / tau
    if verbose:
        return {'value': res, 'avg_arate': np.mean(1 / ARiters), 'ARiters': ARiters}
    else:
        return res

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
