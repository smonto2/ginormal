import numpy as np
from ginormal.utils.rtgin1 import rtgin1

def rtgin(size, alpha, mu, tau, sign, algo='hormann', verbose=False):
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
