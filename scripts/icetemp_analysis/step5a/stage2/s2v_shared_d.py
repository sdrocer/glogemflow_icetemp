"""Variant: observation error shared by the months of one borehole at one point (sensor offset), plus
independent reading noise n0. Refit by pair likelihood at every run (same-borehole pairs carry no
sensor error; different-borehole pairs carry (sigma_i^2 + sigma_j^2)/2), then misfit ratios."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta
import s2m_final_fit as F
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
P, R, w, bf = F.P, F.R, F.w, F.bf
vobs_shared = np.where(bf > 0, P.var_obs.to_numpy(), 0.0)
NAMES = F.NAMES + ['n0']
LO = np.append(F.LO, np.log(1e-6)); HI = np.append(F.HI, np.log(1.0)); ISL = np.append(F.ISL, True)


def unpack(x):
    x = np.clip(x, LO, HI)
    return dict(zip(NAMES, np.where(ISL, np.exp(x), x)))


def fit_run(run):
    d2 = (R[run, P.i] - R[run, P.j]) ** 2
    def f(x):
        p = unpack(x)
        g = F.gamma(p) + vobs_shared + p['n0']
        return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))
    x0 = np.array([np.log(3.0), np.log(3.5), 0.27, np.log(7.0), np.log(0.25), np.log(5.4), np.log(0.45), np.log(0.035), np.log(1e-3)])
    best = None
    for _ in range(4):
        r1 = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(LO, HI)), options=dict(maxiter=5000, ftol=1e-14, gtol=1e-10))
        r2 = minimize(f, r1.x, method='Nelder-Mead', options=dict(maxfev=12000, xatol=1e-7, fatol=1e-9, adaptive=True))
        c = r2 if r2.fun < r1.fun else r1
        if best is None or c.fun < best.fun - 1e-9:
            best = c
        else:
            break
    return run, unpack(best.x), -best.fun


if __name__ == '__main__':
    with mp.get_context('fork').Pool(16) as pool:
        res = pool.map(fit_run, range(32))
    T = pd.DataFrame([dict(run=r, cl=cl, **p) for r, p, cl in res]).set_index('run')
    keep = F.mod['design'].keep.to_numpy()
    print(pd.DataFrame({'design 13': T.loc[13], 'median29': T[keep].median(), 'min29': T[keep].min(), 'max29': T[keep].max()}).round(5).to_string())
    T.to_pickle(W + 'shared_fit_runs.pkl')
