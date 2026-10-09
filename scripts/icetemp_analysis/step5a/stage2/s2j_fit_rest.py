"""Depth-correlation lengths (Matern 3/2 vs exponential) for both parts, and the between-borehole
variance u(z) = u0 exp(-z/zu) (pair likelihood on different-borehole pairs, within part fixed)."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L
import s2i_depthcorr as DC
import s2kern2 as K2
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
pw = DC.pw
P, R, rows = DC.P, DC.R, DC.rows


def pacm(z1, z2, l0, l1, nu):
    a, b = l0 + l1 * z1, l0 + l1 * z2
    L2 = 0.5 * (a * a + b * b)
    return np.sqrt(a * b / L2) * K2.matern(np.abs(z1 - z2) / np.sqrt(L2), nu)


def fit_corr(C, M, fn, nu):
    def f(x):
        l0, l1 = np.exp(x[0]), x[1]
        if l1 < 0:
            return 1e9
        cov, v1, v2 = fn(C, l0, l1, nu)
        t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
        mc = t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())
        return float(np.sum(M['n'] * (M['corr'] - mc) ** 2))
    best = min((minimize(f, [np.log(s0), s1], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxfev=4000))
                for s0 in (0.3, 1, 3, 8) for s1 in (0.0, 0.3, 1.0, 2.0)), key=lambda r: r.fun)
    return np.exp(best.x[0]), best.x[1], best.fun


def within_fn(C, l0, l1, nu):
    v1 = pw['a'] * np.exp(-C.z1 / pw['zs']) * np.where(C.ice1, pw['q'], 1)
    v2 = pw['a'] * np.exp(-C.z2 / pw['zs']) * np.where(C.ice2, pw['q'], 1)
    omr = 1 - (pw['w'] * np.cos(W12 * C.lag) + (1 - pw['w']) * np.exp(-C.lag / pw['tau']))
    return 2 * np.sqrt(v1 * v2) * pacm(C.z1.to_numpy(), C.z2.to_numpy(), l0, l1, nu) * omr, 2 * v1 * omr, 2 * v2 * omr


def between_fn(C, l0, l1, nu):
    return pacm(C.z1.to_numpy(), C.z2.to_numpy(), l0, l1, nu), np.ones(len(C)), np.ones(len(C))


res = {}
for label, same, fn in (('within', True, within_fn), ('between', False, between_fn)):
    C = DC.build(same); M = DC.measured(C)
    for nu in (0.5, 1.5):
        l0, l1, sse = fit_corr(C, M, fn, nu)
        res[(label, nu)] = (l0, l1, sse)
        print(f'{label:8s} nu {nu}: l(z) = {l0:.3f} + {l1:.3f} z   SSE {sse:.1f}')
    l0, l1, _ = res[(label, 1.5)]
    cov, v1, v2 = fn(C, l0, l1, 1.5)
    t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
    M['matern'] = t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())
    print(M.round(3).to_string())
    res[(label, 'table')] = M
# between-borehole variance by pair likelihood on different-borehole pairs (within part at their lag)
Q = P[~P.same_bh].reset_index(drop=True)
bh = rows.bh_label.to_numpy()
Bi = [set(str(bh[i]).split('/')) for i in Q.i]; Bj = [set(str(bh[j]).split('/')) for j in Q.j]
def bdist(a, b):
    keys = sorted(a | b)
    x = np.array([1 / len(a) if k in a else 0 for k in keys]); y = np.array([1 / len(b) if k in b else 0 for k in keys])
    return 0.5 * np.sum((x - y) ** 2)          # gamma factor: 0.5 |b_i - b_j|^2
Q['bfac'] = [bdist(a, b) for a, b in zip(Bi, Bj)]
ser = rows.groupby(['t2']).size()
Q['w'] = 2.0 / Q.t2.map(rows[~rows.index.isin(np.load(W + 'spike_rows_month.npy'))].groupby('t2').size()).to_numpy()
z, lag, ice = Q.depth.to_numpy(), Q.dlag.to_numpy(), Q.ice.to_numpy()
gw = pw['a'] * np.exp(-z / pw['zs']) * np.where(ice, pw['q'], 1) * (1 - (pw['w'] * np.cos(W12 * lag) + (1 - pw['w']) * np.exp(-lag / pw['tau'])))
d2 = (R[13, Q.i] - R[13, Q.j]) ** 2
def nll(x, form):
    u = np.exp(x[0]) * (np.exp(-z / np.exp(x[1])) if form == 'exp' else 1.0)
    g = gw + u * 2 * Q.bfac.to_numpy() + Q.var_obs.to_numpy()      # 2 * bfac = 1 for two single boreholes
    return np.sum(Q.w * (0.5 * np.log(2 * g) + d2 / (4 * g)))
for form in ('const', 'exp'):
    best = min((minimize(nll, [np.log(s), np.log(zz)], args=(form,), method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxfev=4000))
                for s in (0.3, 1, 3) for zz in (5, 20, 100)), key=lambda r: r.fun)
    u0 = np.exp(best.x[0]); zu = np.exp(best.x[1]) if form == 'exp' else np.inf
    res[('u', form)] = (u0, zu, -best.fun)
    print(f'between-borehole u(z) {form}: u0 {u0:.3f} zu {zu:.3g}  CL {-best.fun:.2f}  (pairs {len(Q)}, points {Q.t2.nunique()})')
pickle.dump(res, open(W + 'rest_d13.pkl', 'wb'))
