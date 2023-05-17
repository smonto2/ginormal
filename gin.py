import numpy as np
from scipy import special, stats
from cardano_method import CubicEquation

# Generalized inverse normal (log) density
# Standardized density setting tau = 1
def dgin1(z, a, m, log = True, quasi = False):
    # Impose restrictions on parameters (a > 1)
    if (z == 0) | (a <= 1):
        res = -np.inf
    else:
        am = 0.5*(a - 1)
        lkern = (-a)*np.log(np.abs(z)) - 0.5*(1/z - m)**2
        if quasi:
            res = lkern
        else:
            hyp = special.hyp1f1(am, 0.5, 0.5*m**2)
            lcons = -0.5*m**2 + am*np.log(2) + special.loggamma(am) + np.log(hyp)
            res = lkern - lcons
    if log:
        return res
    else:
        return np.exp(res)

# General density with tau from (0, infinity)
def dgin(z, a, m, t, log = True, quasi = False):
    if (t <= 0):
        res = -np.inf
    else:
        mt = m/t
        if quasi:
            res = dgin1(z*t, a, mt, True, True) + a*np.log(t)
        else:
            res = dgin1(z*t, a, mt, True, False) + np.log(t)
    if log:
        return res
    else:
        return np.exp(res)

# Truncated generalized inverse normal log-densities
def pbdam(a, m):
    am = -a + 1
    return special.pbdv(am, m)[0]

def dtgin(z, a, m, t, sign = True, log = True, quasi = False):
    # Sign = True, means truncate to z > 0; Sign = False is to z < 0
    if (sign):
        t_z = (z <= 0)
        mult = -1
    else:
        t_z = (z >= 0)
        mult = 1

    # Impose restrictions and compute density
    if t_z | (a <= 1) | (t <= 0):
        res = -np.inf
    else:
        mt = m/t
        lkern = dgin1(z*t, a, mt, True, True) + np.log(t)
        if quasi:
            res = lkern
        else:
            lcons = -0.25*mt**2 + special.loggamma(a - 1) + np.log(pbdam(a, mult*mt))
            res = lkern - lcons
    if log:
        return res
    else:
        return np.exp(res)

# CDF of generalized inverse normal at 0:
def F0gin(a, m, t, sign = True):
    # If sign = True, return upper area (1 - F(0))
    # sign = False, return lower area F(0)
    mult = 1 - 2*sign
    mt = m/t
    parhyp = np.log(pbdam(a, mult*mt)) - np.log(special.hyp1f1(0.5*(a-1), 0.5, 0.5*mt**2))
    lcons = 0.25*mt**2 + 0.5*(a-3)*np.log(2) - 0.5*np.log(np.pi) + special.loggamma(0.5*a) + parhyp
    return np.exp(lcons)

# Generator of truncated generalized inverse normal variables
def mshift(z, mode, a, m, shift = True):
    if (shift):
        x = z - mode
    else:
        x = 1
    return x*np.sqrt(dgin1(z, a, m, False, True))

def rtgin1(a, m, sign, algo = True, verbose = False):
    # If sign = True, draw from the positive region Z > 0
    # If sign = False, draw from the negative region Z < 0
    # Compute necessary values for ratio-of-uniforms method
    if (sign):
        mode = (-m + np.sqrt(m**2 + 4*a))/(2*a)
        mult = -1
    else:
        mode = (-m - np.sqrt(m**2 + 4*a))/(2*a)
        mult = 1

    # Draw from minimal bounding rectangle and accept draw
    vmax = mshift(mode, mode, a, m, False)
    if algo:
        # Hörmann and Leydold (2014) using Cardano method
        roots = np.real(CubicEquation([2-a, -m + a*mode, 1 + m*mode, -mode]).answers)
        roots = roots[-mult*roots > 0]
        uvals = np.sort([mshift(roots[j], mode, a, m) for j in range(2)])
    else:
        # Leydold (2001) using the proportionality constant
        lcons = -0.25*m**2 + special.loggamma(a - 1) + np.log(pbdam(a, mult*m))
        vp = np.exp(lcons)/vmax
        uvals = np.array([-vp, vp])
    
    # Acceptance-Rejection algorithm (ratio-of-uniforms)
    test = False
    counter = 0
    while not test:
        u = stats.uniform.rvs(uvals[0], uvals[1] - uvals[0], size = 1)[0]
        v = stats.uniform.rvs(0, vmax, size = 1)[0]
        x = (u/v) + mode
        test = (2*np.log(v) <= dgin1(x, a, m, True, True)) & (-mult*x > 0)
        counter += 1
    if verbose:
        return {'value': x, 'ARiters': counter}
    else:
        return x

def rtgin(size, a, m, t, sign, algo = True, verbose = False):
    # If sign = True, draw from the positive region Z > 0
    # If sign = False, draw from the negative region Z < 0
    
    # If algo = True, use the Minimal Bounding Rectangle from Hörmann and Leydold (2014),
    # if algo = False, use the one from Leydold (2001)
    res = np.zeros(size)
    if verbose: ARiters = np.copy(res)
    mt = m/t
    for i in range(size):
        temp = rtgin1(a, mt, sign, algo, verbose)
        if verbose:
            res[i] = temp['value']/t
            ARiters[i] = temp['ARiters']
        else:
            res[i] = temp/t
    if verbose:
        return {'value': res, 'avg_arate': np.mean(1/ARiters), 'ARiters': ARiters}
    else:
        return res

# Generator of generalized inverse normal variables
def rgin(size, a, m, t, sign = True, algo = True):
    # If sign is true, more samples are expected from the positive region (e.g., mu > 0),
    # i.e., P(>0) >= 50%. Probability calculation of that region is more stable
    # If sign is negative, you expect P(<0) >= 50% (e.g., mu < 0) and that region's probability is used
    
    # If algo = True, use the Minimal Bounding Rectangle from Hörmann and Leydold (2014),
    # if algo = False, use the one from Leydold (2001)
    
    F0r = F0gin(a, m, t, sign)
    res = np.zeros(size)
    side = (stats.uniform.rvs(size = size) <= F0r)
    sum_side = sum(side)
    res[side == 0] = rtgin(size - sum_side, a, m, t, not sign, algo)
    res[side == 1] = rtgin(sum_side, a, m, t, sign, algo)
    return res
