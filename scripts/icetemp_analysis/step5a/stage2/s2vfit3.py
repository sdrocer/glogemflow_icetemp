"""Semivariogram with a knot-wise (log-linear) seasonal variance in depth, firn/ice factor,
periodic + exponential memory, optional slower interannual part and point-centred drift."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
from scipy.optimize import minimize

W12 = 2 * np.pi / 12.0
KNOTS = np.array([1.0, 4.0, 10.0, 20.0, 40.0])


def v_depth(z, logv, knots=KNOTS):
    """Variance at depth z: log-linear between knots, constant beyond the end knots."""
    return np.exp(np.interp(z, knots, logv))


def gamma_model(p, z, lag, ice):
    logv = np.array([p[f'v{k}'] for k in range(len(KNOTS))])
    v = v_depth(z, logv) * np.where(ice, p.get('q', 1.0), 1.0)
    w = p.get('w', 0.0)
    rho = w * np.cos(W12 * lag) + (1 - w) * np.exp(-lag / p['tau'])
    g = v * (1 - rho)
    if 'ay' in p:   # slower interannual part, same depth shape, Gaussian in time
        g = g + p['ay'] * v * (1 - np.exp(-(lag / p['ty']) ** 2))
    if 'beta' in p:
        g = g + 0.5 * (p['beta'] * lag / 12.0) ** 2 * np.where(ice, p.get('qd', 1.0), 1.0)
    return g


SPEC = {'w': (0.0, 1.0, 0), 'tau': (0.05, 600, 1), 'q': (1e-3, 5, 1), 'beta': (1e-5, 1.0, 1),
        'qd': (1e-3, 5, 1), 'ay': (1e-4, 20, 1), 'ty': (6, 2000, 1)}
for k in range(len(KNOTS)):
    SPEC[f'v{k}'] = (np.log(1e-5), np.log(50.0), 0)      # already log variance


def negcl(p, data):
    z, lag, ice, d2, vobs, w = data
    g = gamma_model(p, z, lag, ice) + vobs
    return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))


def fit(start, data, fixed=None):
    fixed = fixed or {}
    names = [k for k in start if k not in fixed]
    isl = np.array([SPEC[k][2] for k in names], bool)
    lo = np.array([np.log(SPEC[k][0]) if l else SPEC[k][0] for k, l in zip(names, isl)])
    hi = np.array([np.log(SPEC[k][1]) if l else SPEC[k][1] for k, l in zip(names, isl)])

    def unpack(x):
        x = np.clip(x, lo, hi)
        return dict(fixed, **dict(zip(names, np.where(isl, np.exp(x), x))))
    x0 = np.clip(np.array([np.log(start[k]) if l else start[k] for k, l in zip(names, isl)]), lo, hi)
    f = lambda x: negcl(unpack(x), data)
    best = None
    for _ in range(4):
        res = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                       options=dict(maxiter=4000, ftol=1e-14, gtol=1e-10))
        res2 = minimize(f, res.x, method='Nelder-Mead', options=dict(maxfev=8000, xatol=1e-7, fatol=1e-9, adaptive=True))
        cand = res2 if res2.fun < res.fun else res
        if best is None or cand.fun < best.fun - 1e-9:
            best = cand
        else:
            break
    return unpack(best.x), -best.fun, len(names)
