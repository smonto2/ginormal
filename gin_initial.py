import numpy as np
from scipy import special, stats
from cardano_method import CubicEquation

# Generalized inverse normal (log) density
# Standardized density setting tau = 1
def dgin1(z, a, m, log=True, quasi=False):
    # Impose restrictions on parameters (a > 1)
    if (z == 0) | (a <= 1):
        res = -np.inf
    else:
        am = 0.5 * (a - 1)
        lkern = (-a) * np.log(np.abs(z)) - 0.5 * (1 / z - m) ** 2
        if quasi:
            res = lkern
        else:
            hyp = special.hyp1f1(am, 0.5, 0.5 * m ** 2)
            lcons = -0.5 * m ** 2 + am * np.log(2) + special.loggamma(am) + np.log(hyp)
            res = lkern - lcons
    if log:
        return res
    else:
        return np.exp(res)

# General density with tau from (0, infinity)
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
        raise ValueError("alpha must be more than 1")
    if t <= 0:
        raise ValueError("tau must be more than 0")
    if (t <= 0):
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

# CDF of generalized inverse normal at 0:
def F0gin(a, m, t, sign=True):
    # If sign = True, return upper area (1 - F(0))
    # sign = False, return lower area F(0)
    mult = 1 - 2 * sign
    mt = m / t
    parhyp = np.log(pbd.pbdam(a, mult * mt)) - np.log(special.hyp1f1(0.5 * (a - 1), 0.5, 0.5 * mt ** 2))
    lcons = 0.25 * mt ** 2 + 0.5 * (a - 3) * np.log(2) - 0.5 * np.log(np.pi) + special.loggamma(0.5 * a) + parhyp
    return np.exp(lcons)

# Generator of truncated generalized inverse normal variables
def mshift(z, mode, a, m, shift=True):
    if (shift):
        x = z - mode
    else:
        x = 1
    return x * np.sqrt(dg_1.dgin1(z, a, m, False, True))

def rtgin1(a, m, sign, algo='hormann', verbose=False):
    # If sign = True, draw from the positive region Z > 0
    # If sign = False, draw from the negative region Z < 0
    # Compute necessary values for ratio-of-uniforms method
    if (sign):
        mode = (-m + np.sqrt(m ** 2 + 4 * a)) / (2 * a)
        mult = -1
    else:
        mode = (-m - np.sqrt(m ** 2 + 4 * a)) / (2 * a)
        mult = 1
    if algo == 'hormann':
        algo_1 = True
    elif algo == 'leydold':
        algo_1 = False
    else:
        raise ValueError("algo_method must be either 'hormann' or 'leydold'")

    # Draw from minimal bounding rectangle and accept draw
    vmax = msh.mshift(mode, mode, a, m, False)
    if algo_1:
        # Hörmann and Leydold (2014) using Cardano method
        roots = np.real(CubicEquation([2 - a, -m + a * mode, 1 + m * mode, -mode]).answers)
        roots = roots[-mult * roots > 0]
        uvals = np.sort([msh.mshift(roots[j], mode, a, m) for j in range(2)])
    else:
        # Leydold (2001) using the proportionality constant
        lcons = -0.25 * m ** 2 + special.loggamma(a - 1) + np.log(pbd.pbdam(a, mult * m))
        vp = np.exp(lcons) / vmax
        uvals = np.array([-vp, vp])

    # Acceptance-Rejection algorithm (ratio-of-uniforms)
    test = False
    counter = 0
    max_iter = 100
    while (not test):
        u = stats.uniform.rvs(uvals[0], uvals[1] - uvals[0], size=1)[0]
        v = stats.uniform.rvs(0, vmax, size=1)[0]
        x = (u / v) + mode
        test = (2 * np.log(v) <= dg_1.dgin1(x, a, m, True, True)) & (-mult * x > 0)
        counter += 1
        if counter > max_iter:
            x = mode
            warnings.warn("No candidate draw found for these parameter values. Giving up an returning the mode of the distribution")
            break
    if verbose:
        return {'value': x, 'ARiters': counter}
    else:
        return x

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
    algo : boolean, optional
        algo = True, use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
        algo = False, use the one from Leydold (2001)
        Defaults to True.
    verbose : boolean, optional
        Should the acceptance rate from the ratio-of-uniforms method be provided
        along with additional information? Defaults to False. 
    """
    if a <= 2:
        raise ValueError("alpha must be more than 2")
    if t <= 0:
        raise ValueError("tau must be more than 0")
    if algo == 'hormann':
        algo_1 = True
    elif algo == 'leydold':
        algo_1 = False
    else:
        raise ValueError("algo_method must be either 'hormann' or 'leydold'")
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
    If more samples are expected from the positive region, probability of positive region is used.
    As m controls asymmetry, when m > 0, P(z>0) >= 50%, and this probability is computed.
    If m < 0, P(z<0) >= 50% and this region's probability is used. 

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
    algo : boolean, optional
        algo = True, use the Minimal Bounding Rectangle from Hörmann and Leydold (2014)
        algo = False, use the one from Leydold (2001)
        Defaults to True.
    """
    if a <= 2:
        raise ValueError("alpha must be more than 2")
    if t <= 0:
        raise ValueError("tau must be more than 0")
    sign = (m >= 0)
    F0r = F0gin(a, m, t, sign)
    res = np.zeros(size)
    side = (stats.uniform.rvs(size=size) <= F0r)
    sum_side = sum(side)
    res[side == 0] = rtgin(size - sum_side, a, m, t, not sign, algo)
    res[side == 1] = rtgin(sum_side, a, m, t, sign, algo)
    return res
