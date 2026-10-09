"""Per-fold refit of the Stage-2 month-to-month term K (pair likelihood of s2m_final_fit.py) with the
held-out glacier system's pairs removed, so the T3 held-out RMSE does not use the left-out system
to set its error model. Same likelihood, weights, bounds and starts as Stage 2; fitted at each of
the 29 runs with f_ref >= 0.1 and summarised by the median, as SEASONAL was. The depth-correlation
lengths (lz, lens, ld) stay at their Stage-2 values. A factor without any pair in the training
systems is unidentified there and set to its neutral value (hf = 1, u0 = 0); the fold's rows that
it would act on are listed. Output: kfold.json (per fold and 'all': medians, ranges, pair counts)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import json
import sys
import time
import multiprocessing as mp

import numpy as np
import pandas as pd

S2 = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
OUT = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/revise/'
sys.path.insert(0, S2)
sys.path.insert(0, '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/scripts/icetemp_analysis')
from scipy.optimize import minimize  # noqa: E402
import target_analysis as ta  # noqa: E402
import s2lib as L  # noqa: E402

W12 = 2 * np.pi / 12
obs = ta.load_obs()
mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(S2 + 'spike_rows_month.npy')
R = Y - y[None, :]
P0 = pd.read_pickle(S2 + 'pairs.pkl')
P_all = P0[~(P0.i.isin(spike) | P0.j.isin(spike))].reset_index(drop=True)
ok = np.ones(len(rows), bool)
ok[spike] = False
pts = rows[ok].groupby('t2').size()
bh = rows.bh_label.to_numpy()


def bdist(a, b):
    a, b = set(str(a).split('/')), set(str(b).split('/'))
    keys = sorted(a | b)
    x = np.array([1 / len(a) if k in a else 0 for k in keys])
    yv = np.array([1 / len(b) if k in b else 0 for k in keys])
    return float(np.sum((x - yv) ** 2))


NAMES = ['a', 'zs', 'w', 'tau', 'q', 'hf', 'u0', 'beta']
LO = np.array([np.log(1e-4), np.log(0.2), 0.0, np.log(0.05), np.log(1e-3), np.log(1e-2), np.log(1e-4), np.log(1e-5)])
HI = np.array([np.log(50), np.log(300), 1.0, np.log(600), np.log(5), np.log(100), np.log(20), np.log(1.0)])
ISL = np.array([1, 1, 0, 1, 1, 1, 1, 1], bool)


def unpack(x):
    x = np.clip(x, LO, HI)
    return dict(zip(NAMES, np.where(ISL, np.exp(x), x)))


def setup(held):
    P = P_all if held == 'all' else P_all[P_all.system_short != held].reset_index(drop=True)
    w = 2.0 / P.t2.map(pts).to_numpy()
    bf = np.array([bdist(bh[i], bh[j]) for i, j in zip(P.i, P.j)])
    ice = rows.firn_label.to_numpy()[P.i] == 'ice'
    high = (rows.band.to_numpy()[P.i] >= 4000) & ~ice
    return dict(P=P, w=w, bf=bf, ice=ice, high=high, z=P.depth.to_numpy(), lag=P.dlag.to_numpy(),
                vobs=P.var_obs.to_numpy())


def gamma(p, s):
    m = np.where(s['ice'], p['q'], np.where(s['high'], p['hf'], 1.0))
    g = m * p['a'] * np.exp(-s['z'] / p['zs']) * (p['w'] * (1 - np.cos(W12 * s['lag']))
                                                  + (1 - p['w']) * (1 - np.exp(-s['lag'] / p['tau'])))
    return g + p['u0'] * s['bf'] / 2 + 0.5 * (p['beta'] * s['lag'] / 12.0) ** 2


SETUPS = {}


def fit_one(task):
    held, run = task
    s = SETUPS[held]
    P = s['P']
    d2 = (R[run, P.i.to_numpy()] - R[run, P.j.to_numpy()]) ** 2

    def f(x):
        g = gamma(unpack(x), s) + s['vobs']
        return np.sum(s['w'] * (0.5 * np.log(2 * g) + d2 / (4 * g)))

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
    return held, run, unpack(best.x), -best.fun


if __name__ == '__main__':
    t0 = time.time()
    systems = sorted(rows.system_short.unique())
    folds = ['all'] + sorted(set(obs['T2'].system_short.unique()))
    for h in folds:
        SETUPS[h] = setup(h)
    keep = mod['design'].keep.to_numpy()
    runs = [int(i) for i in np.flatnonzero(keep)]
    tasks = [(h, r) for h in folds for r in runs]
    with mp.get_context('fork').Pool(12) as pool:
        res = pool.map(fit_one, tasks, chunksize=1)
    T = pd.DataFrame([dict(held=h, run=r, cl=cl, **p) for h, r, p, cl in res])
    T.to_pickle(OUT + 'kfold_runs.pkl')
    out = {}
    for h in folds:
        t = T[T.held == h].set_index('run')[NAMES]
        s = SETUPS[h]
        med = t.median()
        info = dict(pairs=int(len(s['P'])), pairs_high_firn=int(s['high'].sum()), pairs_ice=int(s['ice'].sum()),
                    pairs_two_boreholes=int((s['bf'] > 0).sum()))
        neutral = {}
        if info['pairs_high_firn'] == 0:
            neutral['hf'] = 1.0
        if info['pairs_two_boreholes'] == 0:
            neutral['u0'] = 0.0
        if info['pairs_ice'] == 0:
            neutral['q'] = 1.0
        med = med.to_dict()
        med.update(neutral)
        out[h] = dict(median=med, min=t.min().to_dict(), max=t.max().to_dict(), neutral=list(neutral), **info,
                      seasonal=dict(a=med['a'], zs=med['zs'], w=med['w'], tau=med['tau'], ice=med['q'], high=med['hf'],
                                    u0=med['u0'], beta=med['beta']))
        print(f"{h:16s} pairs {info['pairs']:5d} (high firn {info['pairs_high_firn']}, ice {info['pairs_ice']}, two boreholes "
              f"{info['pairs_two_boreholes']}); neutral {list(neutral)}; median " +
              ' '.join(f'{k} {v:.4g}' for k, v in med.items()), flush=True)
    with open(OUT + 'kfold.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(f'done in {time.time() - t0:.0f} s')
