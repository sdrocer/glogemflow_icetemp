"""Depth length growing exponentially with depth, l(z) = l0 exp(z / zl), fitted to the measured
same-borehole increment correlations (Matern 3/2, final variance and time parameters)."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import s2kern2 as K2
import s2n_final_eval as FE      # reuses C, M, within variance setup (prints its own fit again)
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
C, M, pf = FE.C, FE.M, FE.pf


def fn(l0, zl):
    m1 = np.where(C.ice1, pf['q'], np.where(FE.high1, pf['hf'], 1.0)); m2 = np.where(C.ice2, pf['q'], np.where(FE.high2, pf['hf'], 1.0))
    v1 = m1 * pf['a'] * np.exp(-C.z1 / pf['zs']); v2 = m2 * pf['a'] * np.exp(-C.z2 / pf['zs'])
    omr = 1 - (pf['w'] * np.cos(FE.W12 * C.lag) + (1 - pf['w']) * np.exp(-C.lag / pf['tau']))
    dr = 0.5 * (pf['beta'] * C.lag / 12) ** 2
    a, b = l0 * np.exp(C.z1.to_numpy() / zl), l0 * np.exp(C.z2.to_numpy() / zl)
    L2 = 0.5 * (a * a + b * b)
    c = np.sqrt(a * b / L2) * K2.matern(np.abs(C.z1 - C.z2).to_numpy() / np.sqrt(L2), 1.5)
    cd = K2.matern(np.abs(C.z1 - C.z2).to_numpy() / 50.0, 1.5)
    return 2 * np.sqrt(v1 * v2) * c * omr + 2 * dr * cd, 2 * v1 * omr + 2 * dr, 2 * v2 * omr + 2 * dr


def binned(l0, zl):
    cov, v1, v2 = fn(l0, zl)
    t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
    return t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())


def sse(x):
    return float(np.sum(M['n'] * (M['corr'] - binned(np.exp(x[0]), np.exp(x[1]))) ** 2))


best = min((minimize(sse, [np.log(a), np.log(b)], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxfev=4000))
            for a in (0.3, 0.8, 2.0) for b in (2.0, 5.0, 10.0)), key=lambda r: r.fun)
l0, zl = float(np.exp(best.x[0])), float(np.exp(best.x[1]))
print(f'\nl(z) = {l0:.3f} exp(z / {zl:.2f})  SSE {best.fun:.1f}   [linear form SSE 549.3]')
print('length at 1, 2.5, 5, 9, 15, 20, 30 m:', [round(l0 * np.exp(zz / zl), 2) for zz in (1, 2.5, 5, 9, 15, 20, 30)])
M2 = M[['n', 'corr', 'corr_med']].copy()
M2['model_exp_growth'] = binned(l0, zl)
M2['model_linear'] = M['model']
print(M2.round(3).to_string())
d = pickle.load(open(W + 'final_kernel.pkl', 'rb'))
d['final_growth'] = dict(d['final'], l0=l0, zl=zl, lb0=l0, zlb=zl)
d['corr_table_growth'] = M2
pickle.dump(d, open(W + 'final_kernel.pkl', 'wb'))
