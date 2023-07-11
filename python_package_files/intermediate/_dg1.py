# Generalized inverse normal (log) density
# Standardized density setting tau = 1
import numpy as np
from scipy import special

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

