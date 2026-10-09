"""Kernel variants fitted on same-borehole pairs (spike rows out), then misfit ratios with the
point-centred kernel: A pooled, B pooled + high-firn factor (firn at band >= 4000 m), C fitted
without Monte Rosa."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L, s2kern2 as K2
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(W + 'spike_rows_month.npy')
R = Y - y[None, :]
P0 = pd.read_pickle(W + 'pairs.pkl')
P = P0[~(P0.i.isin(spike) | P0.j.isin(spike))].reset_index(drop=True)
ok = np.ones(len(rows), bool); ok[spike] = False
bh = rows.bh_label.to_numpy()
ser = rows[ok].groupby(['t2', 'bh_label']).size()
pts = rows[ok].groupby('t2').size()
P['w_same'] = [2.0 / ser[(t, b)] for t, b in zip(P.t2, bh[P.i])]
P['w_pt'] = 2.0 / P.t2.map(pts).to_numpy()
P['ice'] = rows.firn_label.to_numpy()[P.i] == 'ice'
P['high'] = (rows.band.to_numpy()[P.i] >= 4000) & ~P.ice


def fit_within(Q, run=13, use_hf=False):
    d2 = (R[run, Q.i] - R[run, Q.j]) ** 2
    z, lag, ice, high, vobs, w = (Q.depth.to_numpy(), Q.dlag.to_numpy(), Q.ice.to_numpy(), Q.high.to_numpy(),
                                  Q.var_obs.to_numpy(), Q.w_same.to_numpy())
    names = ['a', 'zs', 'w', 'tau', 'q'] + (['hf'] if use_hf else [])
    lo = np.array([np.log(1e-4), np.log(0.2), 0, np.log(0.05), np.log(1e-3), np.log(1e-2)])[:len(names)]
    hi = np.array([np.log(50), np.log(300), 1, np.log(600), np.log(5), np.log(100)])[:len(names)]
    isl = np.array([1, 1, 0, 1, 1, 1], bool)[:len(names)]

    def unpack(x):
        x = np.clip(x, lo, hi)
        return dict(zip(names, np.where(isl, np.exp(x), x)))

    def f(x):
        p = unpack(x)
        m = np.where(ice, p['q'], 1.0)
        if use_hf:
            m = np.where(high, p['hf'], m)
        g = m * p['a'] * np.exp(-z / p['zs']) * (1 - (p['w'] * np.cos(W12 * lag) + (1 - p['w']) * np.exp(-lag / p['tau']))) + vobs
        return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))
    x0 = np.array([np.log(1.5), np.log(5), 0.3, np.log(4), np.log(0.2), np.log(3)])[:len(names)]
    best = None
    for _ in range(4):
        r1 = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)), options=dict(maxiter=4000, ftol=1e-14, gtol=1e-10))
        r2 = minimize(f, r1.x, method='Nelder-Mead', options=dict(maxfev=8000, xatol=1e-7, fatol=1e-9, adaptive=True))
        c = r2 if r2.fun < r1.fun else r1
        if best is None or c.fun < best.fun - 1e-9:
            best = c
        else:
            break
    p = unpack(best.x)
    if not use_hf:
        p['hf'] = 1.0
    return p, -best.fun


def fit_u(p, run=13):
    Q = P[~P.same_bh].reset_index(drop=True)
    bi = [set(str(bh[i]).split('/')) for i in Q.i]; bj = [set(str(bh[j]).split('/')) for j in Q.j]
    def bfac(a, b):
        keys = sorted(a | b)
        x = np.array([1 / len(a) if k in a else 0 for k in keys]); yv = np.array([1 / len(b) if k in b else 0 for k in keys])
        return np.sum((x - yv) ** 2)            # |b_i - b_j|^2; gamma contribution u * |..|^2 / 2
    bf = np.array([bfac(a, b) for a, b in zip(bi, bj)])
    z, lag, ice, high = Q.depth.to_numpy(), Q.dlag.to_numpy(), Q.ice.to_numpy(), Q.high.to_numpy()
    m = np.where(ice, p['q'], 1.0); m = np.where(high, p['hf'], m)
    gw = m * p['a'] * np.exp(-z / p['zs']) * (1 - (p['w'] * np.cos(W12 * lag) + (1 - p['w']) * np.exp(-lag / p['tau'])))
    d2 = (R[run, Q.i] - R[run, Q.j]) ** 2
    w = Q.w_pt.to_numpy()
    f = lambda x: np.sum(w * (0.5 * np.log(2 * (gw + np.exp(x[0]) * bf / 2 + Q.var_obs)) + d2 / (4 * (gw + np.exp(x[0]) * bf / 2 + Q.var_obs))))
    r = minimize(f, [0.0], method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-10))
    return float(np.exp(r.x[0]))


rest = pickle.load(open(W + 'rest_d13.pkl', 'rb'))
DEPTH = dict(nu=1.5, l0=rest[('within', 1.5)][0], l1=rest[('within', 1.5)][1], lb0=0.05, lb1=rest[('between', 1.5)][1], nub=1.5, centre=True)
same = P[P.same_bh].reset_index(drop=True)
variants = {}
pA, clA = fit_within(same); variants['A pooled'] = dict(pA, u0=fit_u(pA), **DEPTH)
pB, clB = fit_within(same, use_hf=True); variants['B pooled + high-firn factor'] = dict(pB, u0=fit_u(pB), **DEPTH)
pC, clC = fit_within(same[same.system_short != 'Monte Rosa'].reset_index(drop=True)); variants['C without Monte Rosa'] = dict(pC, u0=fit_u(pC), **DEPTH)
print(f'pair likelihood A {clA:.2f} (6 params incl. none), B {clB:.2f} (+hf)')
for k, v in variants.items():
    print(f'{k:30s}', {a: round(float(b), 4) for a, b in v.items() if a not in ('nu', 'nub', 'centre', 'lb0')})
pickle.dump(variants, open(W + 'variants.pkl', 'wb'))


def evaluate(p, run=13, names=('T2', 'T3', 'T3r', 'T3m')):
    res = {}
    for name in names:
        t = ta.target('T3' if name == 'T3m' else name, obs, mod)
        rr, Yt, yt = t.rows, t.Y, t.y
        if name == 'T3m':
            m = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[m].reset_index(drop=True), Yt[:, m], yt[m]
        Sig = ta.error_cov(rr) + K2.seasonal_matrix(rr, p, obs=obs)
        if name == 'T3r':
            Sig = Sig + np.diag(K2.date_spread(rr, obs, mod['Y3'][run]))
        mm = K2.misfit(Yt[run] - yt, Sig, rr)
        pe = mm['per_entity']
        pe['system'] = rr.groupby('e').system_short.first().reindex(pe.e).to_numpy()
        g = pe.groupby('system').agg(n=('n', 'sum'), p=('p', 'sum'), chi2=('chi2', 'sum'), level=('level', 'sum'))
        g['level_ratio'] = g.level / g.p
        g['contrast_ratio'] = (g.chi2 - g.level) / (g.n - g.p).clip(lower=1)
        g['contrast_dims'] = g.n - g.p
        res[name] = (mm, g)
    return res


for k, p in variants.items():
    res = evaluate(p)
    print(f'\n=== {k} (design 13) ===')
    for name, (mm, g) in res.items():
        print(f"{name:4s} n {mm['n']:5d} ratio {mm['ratio']:.3f} level {mm['level']:.3f} contrast {mm['contrast']:.3f}")
    print('T3r by system:'); print(res['T3r'][1][['n', 'p', 'contrast_dims', 'level_ratio', 'contrast_ratio']].round(3).to_string())
