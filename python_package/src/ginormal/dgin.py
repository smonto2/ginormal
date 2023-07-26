import numpy as np
import warnings
from ginormal.utils.dgin1 import dgin1

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
