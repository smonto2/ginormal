import numpy as np
from intermediate import _dg1 as dg_1

# Generator of truncated generalized inverse normal variables
def mshift(z, mode, a, m, shift=True):
    if (shift):
        x = z - mode
    else:
        x = 1
    return x * np.sqrt(dg_1.dgin1(z, a, m, False, True))
