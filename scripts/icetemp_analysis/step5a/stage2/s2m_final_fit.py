"""Joint pair-likelihood fit of the month-to-month kernel on all same-point pairs of T3 month rows
(spike rows out; weight 2 / rows of the point):
  gamma = m(z) a exp(-z/zs) [w (1 - cos(2 pi dt/12)) + (1 - w)(1 - exp(-dt/tau))]       within one borehole
        + u0 |b_i - b_j|^2 / 2                                                              between boreholes
        + (beta dt / 12)^2 / 2                                                              slow drift
  m = q for ice, hf for firn at band >= 4000 m, else 1. Fitted at every run (median reported)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(W + 'spike_rows_month.npy')
R = Y - y[None, :]
P0 = pd.read_pickle(W + 'pairs.pkl')
P = P0[~(P0.i.isin(spike) | P0.j.isin(spike))].reset_index(drop=True)
ok = np.ones(len(rows), bool); ok[spike] = False
pts = rows[ok].groupby('t2').size()
w = 2.0 / P.t2.map(pts).to_numpy()
bh = rows.bh_label.to_numpy()


def bdist(a, b):
    a, b = set(str(a).split('/')), set(str(b).split('/'))
    keys = sorted(a | b)
    x = np.array([1 / len(a) if k in a else 0 for k in keys]); yv = np.array([1 / len(b) if k in b else 0 for k in keys])
    return float(np.sum((x - yv) ** 2))


bf = np.array([bdist(bh[i], bh[j]) for i, j in zip(P.i, P.j)])
ice = rows.firn_label.to_numpy()[P.i] == 'ice'
high = (rows.band.to_numpy()[P.i] >= 4000) & ~ice
z, lag, vobs = P.depth.to_numpy(), P.dlag.to_numpy(), P.var_obs.to_numpy()
NAMES = ['a', 'zs', 'w', 'tau', 'q', 'hf', 'u0', 'beta']
LO = np.array([np.log(1e-4), np.log(0.2), 0.0, np.log(0.05), np.log(1e-3), np.log(1e-2), np.log(1e-4), np.log(1e-5)])
HI = np.array([np.log(50), np.log(300), 1.0, np.log(600), np.log(5), np.log(100), np.log(20), np.log(1.0)])
ISL = np.array([1, 1, 0, 1, 1, 1, 1, 1], bool)


def unpack(x, names=NAMES):
    x = np.clip(x, LO, HI)
    return dict(zip(names, np.where(ISL, np.exp(x), x)))


def gamma(p, z=z, lag=lag, ice=ice, high=high, bf=bf):
    m = np.where(ice, p['q'], np.where(high, p['hf'], 1.0))
    g = m * p['a'] * np.exp(-z / p['zs']) * (p['w'] * (1 - np.cos(W12 * lag)) + (1 - p['w']) * (1 - np.exp(-lag / p['tau'])))
    return g + p['u0'] * bf / 2 + 0.5 * (p['beta'] * lag / 12.0) ** 2


def fit_run(run):
    d2 = (R[run, P.i] - R[run, P.j]) ** 2
    f = lambda x: np.sum(w * (0.5 * np.log(2 * (gamma(unpack(x)) + vobs)) + d2 / (4 * (gamma(unpack(x)) + vobs))))
    x0 = np.array([np.log(3.0), np.log(3.5), 0.27, np.log(9.0), np.log(0.27), np.log(6.0), np.log(0.6), np.log(0.03)])
    best = None
    for _ in range(4):
        r1 = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(LO, HI)),
                      options=dict(maxiter=5000, ftol=1e-14, gtol=1e-10))
        r2 = minimize(f, r1.x, method='Nelder-Mead', options=dict(maxfev=10000, xatol=1e-7, fatol=1e-9, adaptive=True))
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
    keep = mod['design'].keep.to_numpy()
    print(T.loc[13].round(4).to_string())
    print('\nmedian over the 29 runs with f_ref >= 0.1, and range:')
    print(pd.DataFrame({'median29': T[keep].median(), 'min29': T[keep].min(), 'max29': T[keep].max(), 'median32': T.median()}).round(4).to_string())
    T.to_pickle(W + 'final_fit_runs.pkl')
