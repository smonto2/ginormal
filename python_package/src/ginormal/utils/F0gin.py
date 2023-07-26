import numpy as np
from scipy import special
from ginormal.utils.pbdam import pbdam

# CDF of generalized inverse normal at 0:
def F0gin(a, m, t, sign=True):
    # If sign = True, return upper area (1 - F(0))
    # sign = False, return lower area F(0)
    mult = 1 - 2 * sign
    mt = m / t
    parhyp = np.log(pbdam(a, mult * mt)) - np.log(special.hyp1f1(0.5 * (a - 1), 0.5, 0.5 * mt ** 2))
    lcons = 0.25 * mt ** 2 + 0.5 * (a - 3) * np.log(2) - 0.5 * np.log(np.pi) + special.loggamma(0.5 * a) + parhyp
    return np.exp(lcons)


