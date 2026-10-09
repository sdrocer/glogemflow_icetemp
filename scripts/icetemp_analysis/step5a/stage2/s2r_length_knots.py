"""Depth length from knots (2.5, 9, 20, 40 m; log-linear), fitted to the measured same-borehole
increment correlations; bin weights min(n, 300) so deep bins count. Then T3r contrast by depth."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import s2kern2 as K2
import s2n_final_eval as FE
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
C, M, pf = FE.C, FE.M, FE.pf
wt = np.minimum(M['n'].to_numpy(), 300)


def fn(lk):
    m1 = np.where(C.ice1, pf['q'], np.where(FE.high1, pf['hf'], 1.0)); m2 = np.where(C.ice2, pf['q'], np.where(FE.high2, pf['hf'], 1.0))
    v1 = m1 * pf['a'] * np.exp(-C.z1 / pf['zs']); v2 = m2 * pf['a'] * np.exp(-C.z2 / pf['zs'])
    omr = 1 - (pf['w'] * np.cos(FE.W12 * C.lag) + (1 - pf['w']) * np.exp(-C.lag / pf['tau']))
    dr = 0.5 * (pf['beta'] * C.lag / 12) ** 2
    a = K2.length_z(C.z1.to_numpy(), 0, lk=lk); b = K2.length_z(C.z2.to_numpy(), 0, lk=lk)
    L2 = 0.5 * (a * a + b * b)
    c = np.sqrt(a * b / L2) * K2.matern(np.abs(C.z1 - C.z2).to_numpy() / np.sqrt(L2), 1.5)
    cd = K2.matern(np.abs(C.z1 - C.z2).to_numpy() / 50.0, 1.5)
    return 2 * np.sqrt(v1 * v2) * c * omr + 2 * dr * cd, 2 * v1 * omr + 2 * dr, 2 * v2 * omr + 2 * dr


def binned(lk):
    cov, v1, v2 = fn(lk)
    t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
    return t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())


def sse(x, w=wt):
    return float(np.sum(w * (M['corr'].to_numpy() - binned(np.exp(x)).to_numpy()) ** 2))


best = min((minimize(sse, np.log(s0), method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxfev=8000))
            for s0 in ([1.2, 4, 15, 50], [2, 6, 40, 100], [1, 3, 10, 30])), key=lambda r: r.fun)
lk = np.exp(best.x)
print('knot lengths (2.5, 9, 20, 40 m):', lk.round(2), f'capped-weight SSE {best.fun:.1f}',
      f'(linear form, same weights: {sse(np.log(K2.length_z(K2.LKNOTS, FE.FINAL["l0"], FE.FINAL["l1"]))):.1f})')
M3 = M[['n', 'corr', 'corr_med']].copy()
M3['model_knots'] = binned(lk).to_numpy()
M3['model_linear'] = M['model']
print(M3.round(3).to_string())
d = pickle.load(open(W + 'final_kernel.pkl', 'rb'))
d['final_knots'] = dict(d['final'], lk=list(lk), lkb=list(lk))
d['corr_table_knots'] = M3
pickle.dump(d, open(W + 'final_kernel.pkl', 'wb'))
