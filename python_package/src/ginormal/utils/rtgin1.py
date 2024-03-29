import numpy as np
from scipy import special, stats
from cardano_method import CubicEquation
import warnings

from ginormal.utils.mshift import mshift
from ginormal.utils.dgin1 import dgin1
from ginormal.utils.pbdam import pbdam

def rtgin1(a, m, sign, algo='hormann', verbose=False):
    # If sign = True, draw from the positive region Z > 0
    # If sign = False, draw from the negative region Z < 0
    if algo == 'hormann':
        algo_1 = True
    elif algo == 'leydold':
        algo_1 = False
    else:
        raise ValueError("algo_method must be either 'hormann' or 'leydold'")
    
    # Compute necessary values for ratio-of-uniforms method
    mult = 1 - 2*sign
    mode = (-m - mult*np.sqrt(m ** 2 + 4 * a)) / (2 * a)
    
    # Draw from minimal bounding rectangle and accept draw
    vmax = mshift(mode, mode, a, m, False)
    if algo_1:
        # Hörmann and Leydold (2014) using Cardano method
        roots = np.real(CubicEquation([2 - a, -m + a * mode, 1 + m * mode, -mode]).answers)
        roots = roots[-mult * roots > 0]
        uvals = np.sort([mshift(roots[j], mode, a, m) for j in range(2)])
    else:
        # Leydold (2001) using the proportionality constant
        lcons = -0.25 * m ** 2 + special.loggamma(a - 1) + np.log(pbdam(a, mult * m))
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
        test = (2 * np.log(v) <= dgin1(x, a, m, True, True)) & (-mult * x > 0)
        counter += 1
        if counter > max_iter:
            x = mode
            warnings.warn("No candidate draw found for these parameter values. Giving up an returning the mode of the distribution")
            break
    if verbose:
        return {'value': x, 'ARiters': counter}
    else:
        return x
