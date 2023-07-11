import numpy as np
import warnings
from scipy import special, stats
from intermediate._dg1 import dgin1 as dgin1
from intermediate._rtg1 import rtgin1 as rtgin1
from intermediate._F0g import F0gin as F0gin
from intermediate._pbd import pbdam as pbdam

# General density with tau from (0, infinity)
def dgin(z, a, m, t, log=True, quasi=False):
    """
    Density for the generalized inverse normal distribution.

    Parameters
    ----------
    z : real number
        Distribution is supported on real line, z can be any real number.
        Currently, only scalar z is supported. Density is set to be 0 at z = 0.
    a : real number > 1
        Degrees of freedom parameter. Function returns 0 density if a <= 1.
    m : real number
        Similar to location parameter, controls asymmetry of the distribution.
    t : real number > 0
        Similar to scale parameter, controls spread of the distribution.
    log : boolean, optional
        Should the log of the density be returned? Default is True.
    quasi : boolean, optional
        The quasi-density or kernel is the density without the normalization constant.
        Should the quasi-density value be returned? Default is False.
    """
    if a <= 1:
        warnings.warn("a should be greater than 1")
        res = -np.inf
    elif t <= 0:
        warnings.warn("t should be greater than 0")
        res = -np.inf
    elif z == 0:
        res = -np.inf
    else:
        mt = m / t
        if quasi:
            res = dgin1(z * t, a, mt, True, True) + a * np.log(t)
        else:
            res = dgin1(z * t, a, mt, True, False) + np.log(t)
    
    if log:
        return res
    else:
        return np.exp(res)

def dtgin(z, a, m, t, sign=True, log=True, quasi=False):    
    """
    Density for the generalized inverse normal distribution truncated to
    the positive or negative reals.

    Parameters
    ----------
    z : real number
        Distribution is supported on real line, z can be any real number.
        Currently, only scalar z is supported. Density is set to be 0 at z = 0.
    a : real number > 1
        Degrees of freedom parameter. Function returns 0 density if a <= 1.
    m : real number
        Similar to location parameter, controls asymmetry of the distribution.
    t : real number > 0
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
    if a <= 1:
        warnings.warn("a should be greater than 1")
        res = -np.inf
    elif t <= 0:
        warnings.warn("t should be greater than 0")
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
        mt = m / t
        lkern = dgin1(z * t, a, mt, True, True) + np.log(t)
        if quasi:
            res = lkern
        else:
            mult = 1 - 2*sign
            lcons = -0.25 * mt ** 2 + special.loggamma(a - 1) + np.log(pbdam(a, mult * mt))
            res = lkern - lcons
    
    if log:
        return res
    else:
        return np.exp(res)

def rtgin(size, a, m, t, sign=True, algo='hormann', verbose=False):
    """
    Generating random numbers from the generalized inverse normal distribution
    truncated to the positive or negative reals. Currently, only values of a > 2 are supported.
    For Bayesian posterior sampling, a is always larger than 2 even for non-informative priors.

    Parameters
    ----------
    size : integer
        Number of desired draws. Output is numpy vector of length equal to size.
    a : real number > 2
        Degrees of freedom parameter. Currently, only values of a > 2 are supported.
    m : real number
        Similar to location parameter, controls asymmetry of the distribution.
    t : real number > 0
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
    if a <= 2:
        raise ValueError("a should be greater than 2")
    if t <= 0:
        raise ValueError("t should be greater than 0")
    if (algo != 'hormann') & (algo != 'leydold'):
        raise ValueError("algo should be either 'hormann' or 'leydold'")
        
    res = np.zeros(size)
    if verbose: ARiters = np.copy(res)
    mt = m / t
    for i in range(size):
        temp = rtgin1(a, mt, sign, algo, verbose)
        if verbose:
            res[i] = temp['value'] / t
            ARiters[i] = temp['ARiters']
        else:
            res[i] = temp / t
    if verbose:
        return {'value': res, 'avg_arate': np.mean(1 / ARiters), 'ARiters': ARiters}
    else:
        return res


# Generator of generalized inverse normal variables
def rgin(size, a, m, t, algo='hormann'):
    """
    Generating random numbers from the generalized inverse normal distribution
    supported on the real numbers. Currently, only values of a > 2 are supported.
    For Bayesian posterior sampling, a is always larger than 2 even for non-informative priors.
    As m controls asymmetry, when m > 0, P(truncation region) = P(z>0) >= 50%,
    and this probability is computed. If m < 0, P(z<0) >= 50% and this region's probability is used. 

    Parameters
    ----------
    size : integer
        Number of desired draws. Output is numpy vector of length equal to size.
    a : real number > 2
        Degrees of freedom parameter. Currently, only values of a > 2 are supported.
    m : real number
        Similar to location parameter, controls asymmetry of the distribution.
    t : real number > 0
        Similar to scale parameter, controls spread of the distribution.
    algo : string, optional
        algo = 'hormann', use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
        algo = 'leydold', use the one from Leydold (2001)
        Defaults to 'hormann'.
    """
    if a <= 2:
        raise ValueError("a should be greater than 2")
    if t <= 0:
        raise ValueError("t should be greater than 0")
    if (algo != 'hormann') & (algo != 'leydold'):
        raise ValueError("algo should be either 'hormann' or 'leydold'")
    sign = (m >= 0)
    F0r = F0gin(a, m, t, sign)
    res = np.zeros(size)
    side = (stats.uniform.rvs(size=size) <= F0r)
    sum_side = sum(side)
    res[side == 0] = rtgin(size - sum_side, a, m, t, not sign, algo)
    res[side == 1] = rtgin(sum_side, a, m, t, sign, algo)
    return res
