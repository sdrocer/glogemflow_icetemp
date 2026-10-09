"""Fit parametric semivariograms to same-point pairs of T3 month rows (pairwise composite likelihood:
each pair's residual difference ~ N(0, 2 gamma(z, lag)), weight 2 / rows of its point)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, pandas as pd
from scipy.optimize import minimize

W12 = 2 * np.pi / 12.0
# parameter -> (lower, upper, log?)
SPEC = {'ap': (1e-4, 50, 1), 'ar': (1e-4, 50, 1), 'zs': (0.2, 300, 1), 'tr': (0.02, 120, 1), 'tp': (6, 5000, 1),
        'al': (1e-4, 50, 1), 'zl': (0.2, 1000, 1), 'bl': (1e-5, 20, 1), 'tl': (2, 5000, 1),
        'ad': (1e-4, 50, 1), 'td': (12, 5000, 1), 'zd': (0.2, 1000, 1)}


def gamma_model(p, z, lag):
    """Semivariogram between two rows of one point at depth z [m], lag [months] (no obs error)."""
    vs = np.exp(-z / p['zs'])
    g = p.get('ap', 0) * vs * (1 - np.cos(W12 * lag) * (np.exp(-lag / p['tp']) if 'tp' in p else 1.0))
    g = g + p.get('ar', 0) * vs * (1 - np.exp(-lag / p['tr']))
    if 'al' in p:                                     # slow part (gaussian in time)
        vl = p['al'] * np.exp(-z / p['zl']) + p.get('bl', 0.0)
        g = g + vl * (1 - np.exp(-(lag / p['tl']) ** 2))
    if 'ad' in p:                                     # second, decadal slow part
        vd = p['ad'] * (np.exp(-z / p['zd']) if 'zd' in p else 1.0)
        g = g + vd * (1 - np.exp(-(lag / p['td']) ** 2))
    return g


def negcl(p, z, lag, d2, vobs, w):
    g = gamma_model(p, z, lag) + vobs
    return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))


def fit(start, z, lag, d2, vobs, w, fixed=None):
    fixed = fixed or {}
    names = [k for k in start if k not in fixed]
    lo = np.array([np.log(SPEC[k][0]) if SPEC[k][2] else SPEC[k][0] for k in names])
    hi = np.array([np.log(SPEC[k][1]) if SPEC[k][2] else SPEC[k][1] for k in names])
    isl = np.array([SPEC[k][2] for k in names], bool)

    def unpack(x):
        x = np.clip(x, lo, hi)
        return dict(fixed, **dict(zip(names, np.where(isl, np.exp(x), x))))
    x0 = np.clip(np.array([np.log(start[k]) if l else start[k] for k, l in zip(names, isl)]), lo, hi)
    f = lambda x: negcl(unpack(x), z, lag, d2, vobs, w)
    best = None
    for _ in range(3):
        res = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                       options=dict(maxiter=2000, ftol=1e-13, gtol=1e-9))
        res2 = minimize(f, res.x, method='Nelder-Mead', options=dict(maxfev=4000, xatol=1e-6, fatol=1e-8, adaptive=True))
        cand = res2 if res2.fun < res.fun else res
        if best is None or cand.fun < best.fun - 1e-9:
            best = cand
        else:
            break
    return unpack(best.x), -best.fun, len(names)
