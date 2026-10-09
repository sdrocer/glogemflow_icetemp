"""REML fits of the month-to-month (seasonal) error kernel on T3 month rows, one fixed mean per
pooled point (so the level term S, constant within a point, drops out exactly)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import sys; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

W12 = 2 * np.pi / 12.0


def paciorek(z, l0, l1):
    """Non-stationary exponential depth correlation (Paciorek & Schervish 2006), length l0 + l1 z."""
    l = l0 + l1 * z
    L2 = 0.5 * (l[:, None] ** 2 + l[None, :] ** 2)
    return np.sqrt(np.outer(l, l) / L2) * np.exp(-np.abs(z[:, None] - z[None, :]) / np.sqrt(L2))


def kernel_parts(p, b):
    """Month-to-month kernel of one block (precomputed b['dt'], b['dz'], b['zz'], b['cos'], b['z'])."""
    dt, dz, zz = b['dt'], b['dz'], b['zz']
    c = paciorek(b['z'], p['l0'], p.get('l1', 0.0))
    vs = p['a'] * np.exp(-zz / p['ze'])
    w = p.get('w', 0.0)
    per = b['cos'] * (np.exp(-dt / p['tp']) if 'tp' in p else 1.0)
    K = vs * c * (w * per + (1 - w) * np.exp(-dt / p['tr']))
    if 'aC' in p:
        vC = p['aC'] * (np.exp(-zz / p['zC']) if 'zC' in p else 1.0)
        K = K + vC * np.exp(-dz / p.get('lC', 50.0)) * np.exp(-(dt / p['tC']) ** 2)
    return K


class Blocks:
    def __init__(self, rows, sel=None, keep_single=False):
        sel = np.ones(len(rows), bool) if sel is None else np.asarray(sel, bool)
        self.blocks = []
        sub_all = rows[sel]
        pos = np.flatnonzero(sel)
        for e, ii in sub_all.groupby('e').indices.items():
            idx = pos[ii]
            sub = rows.iloc[idx]
            pts, inv = np.unique(sub.t2.to_numpy(), return_inverse=True)
            if len(pts) == len(idx) and not keep_single:
                continue
            Z = np.zeros((len(idx), len(pts))); Z[np.arange(len(idx)), inv] = 1.0
            z = sub.depth.to_numpy(float)
            t = sub.tmon.to_numpy(float) + sub.frac.to_numpy(float)
            dt = np.abs(t[:, None] - t[None, :])
            self.blocks.append(dict(idx=idx, z=z, t=t, s2=sub.sigma.to_numpy(float) ** 2, Z=Z, dt=dt,
                                    dz=np.abs(z[:, None] - z[None, :]), zz=0.5 * (z[:, None] + z[None, :]),
                                    cos=np.cos(W12 * dt)))
        self.n = sum(len(b['idx']) for b in self.blocks)
        self.p = sum(b['Z'].shape[1] for b in self.blocks)


def reml(p, B, r, parts=False):
    tot, quad = 0.0, 0.0
    for b in B.blocks:
        S = kernel_parts(p, b)
        S[np.diag_indices_from(S)] += b['s2']
        try:
            c = cho_factor(S, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return (-np.inf, np.nan) if parts else -np.inf
        rr = r[b['idx']]
        Sir, SiZ = cho_solve(c, rr, check_finite=False), cho_solve(c, b['Z'], check_finite=False)
        A = b['Z'].T @ SiZ
        ca = cho_factor(A, lower=True, check_finite=False)
        u = b['Z'].T @ Sir
        q = rr @ Sir - u @ cho_solve(ca, u, check_finite=False)
        tot += 2 * np.log(np.diag(c[0])).sum() + 2 * np.log(np.diag(ca[0])).sum() + q
        quad += q
    return (-0.5 * tot, quad) if parts else -0.5 * tot


SPEC = {  # name -> (lower, upper, transform)
    'a': (1e-3, 50, 'log'), 'ze': (0.3, 300, 'log'), 'l0': (0.05, 200, 'log'), 'l1': (0.0, 5, 'lin'),
    'w': (0.0, 1.0, 'lin'), 'tp': (1, 5000, 'log'), 'tr': (0.05, 500, 'log'),
    'aC': (1e-4, 20, 'log'), 'zC': (0.3, 1000, 'log'), 'tC': (3, 5000, 'log'), 'lC': (0.5, 1000, 'log'),
}


def fit(B, r, start, fixed=None, maxiter=400):
    fixed = fixed or {}
    names = [k for k in start if k not in fixed]
    islog = np.array([SPEC[k][2] == 'log' for k in names])
    lo = np.array([np.log(SPEC[k][0]) if l else SPEC[k][0] for k, l in zip(names, islog)])
    hi = np.array([np.log(SPEC[k][1]) if l else SPEC[k][1] for k, l in zip(names, islog)])

    def unpack(x):
        x = np.clip(x, lo, hi)
        v = np.where(islog, np.exp(x), x)
        return dict(fixed, **dict(zip(names, v)))

    x0 = np.clip(np.array([np.log(start[k]) if l else start[k] for k, l in zip(names, islog)]), lo, hi)

    def f(x):
        v = -reml(unpack(x), B, r)
        return v if np.isfinite(v) else 1e12

    res = minimize(f, x0, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                   options=dict(maxiter=maxiter, eps=1e-5, ftol=1e-10, gtol=1e-6))
    res2 = minimize(f, res.x, method='Nelder-Mead', options=dict(maxfev=600, xatol=1e-4, fatol=1e-5, adaptive=True))
    best = res2 if res2.fun < res.fun else res
    best2 = minimize(f, best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                     options=dict(maxiter=maxiter, eps=1e-6, ftol=1e-12, gtol=1e-7))
    best = best2 if best2.fun < best.fun else best
    return unpack(best.x), -best.fun, len(names)
