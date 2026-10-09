"""Semivariogram model with firn/ice factor, seasonal (periodic + short noise), interannual and a
point-centred drift; pairwise composite likelihood on same-point pairs."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
from scipy.optimize import minimize

W12 = 2 * np.pi / 12.0
SPEC = {'a': (1e-4, 50, 1), 'zs': (0.2, 300, 1), 'w': (0.0, 1.0, 0), 'tr': (0.02, 120, 1),
        'ay': (1e-5, 50, 1), 'zy': (0.2, 1000, 1), 'ty': (3, 2000, 1),
        'q': (1e-3, 3, 1), 'beta': (1e-5, 1.0, 1), 'qd': (1e-3, 3, 1)}


def gamma_model(p, z, lag, ice):
    """Same-point semivariogram [C^2] at depth z [m], lag [months]; ice: bool per pair."""
    m = np.where(ice, p.get('q', 1.0), 1.0)
    vs = p['a'] * np.exp(-z / p['zs'])
    w = p.get('w', 0.0)
    g = vs * (w * (1 - np.cos(W12 * lag)) + (1 - w) * (1 - np.exp(-lag / p['tr'])))
    if 'ay' in p:
        g = g + p['ay'] * np.exp(-z / p.get('zy', p['zs'])) * (1 - np.exp(-(lag / p['ty']) ** 2))
    g = m * g
    if 'beta' in p:
        g = g + np.where(ice, p.get('qd', 1.0), 1.0) * 0.5 * (p['beta'] * lag / 12.0) ** 2
    return g


def negcl(p, z, lag, ice, d2, vobs, w):
    g = gamma_model(p, z, lag, ice) + vobs
    return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))


def fit(start, data, fixed=None):
    z, lag, ice, d2, vobs, w = data
    fixed = fixed or {}
    names = [k for k in start if k not in fixed]
    isl = np.array([SPEC[k][2] for k in names], bool)
    lo = np.array([np.log(SPEC[k][0]) if l else SPEC[k][0] for k, l in zip(names, isl)])
    hi = np.array([np.log(SPEC[k][1]) if l else SPEC[k][1] for k, l in zip(names, isl)])

    def unpack(x):
        x = np.clip(x, lo, hi)
        return dict(fixed, **dict(zip(names, np.where(isl, np.exp(x), x))))
    x0 = np.clip(np.array([np.log(start[k]) if l else start[k] for k, l in zip(names, isl)]), lo, hi)
    f = lambda x: negcl(unpack(x), z, lag, ice, d2, vobs, w)
    best = None
    for _ in range(4):
        res = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                       options=dict(maxiter=3000, ftol=1e-14, gtol=1e-10))
        res2 = minimize(f, res.x, method='Nelder-Mead', options=dict(maxfev=6000, xatol=1e-7, fatol=1e-9, adaptive=True))
        cand = res2 if res2.fun < res.fun else res
        if best is None or cand.fun < best.fun - 1e-9:
            best = cand
        else:
            break
    return unpack(best.x), -best.fun, len(names)
