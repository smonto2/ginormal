import numpy as np
import warnings
from scipy import special
from ginormal.utils.dgin1 import dgin1
from ginormal.utils.pbdam import pbdam

def dtgin(z, alpha, mu, tau, sign, log=True, quasi=False):    
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
