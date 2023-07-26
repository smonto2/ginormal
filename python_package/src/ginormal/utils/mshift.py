import numpy as np
from ginormal.utils.dgin1 import dgin1

# Generator of truncated generalized inverse normal variables
def mshift(z, mode, a, m, shift=True):
    if (shift):
        x = z - mode
    else:
        x = 1
    return x * np.sqrt(dgin1(z, a, m, False, True))
