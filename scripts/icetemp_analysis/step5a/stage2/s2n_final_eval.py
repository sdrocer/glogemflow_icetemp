"""Final kernel: depth length refit (Matern 3/2, within-borehole increments), then misfit ratios."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L, s2kern2 as K2
import s2i_depthcorr as DC
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
obs, mod = DC.obs, DC.mod
T = pd.read_pickle(W + 'final_fit_runs.pkl')
pf = T.loc[13].drop('cl').to_dict()
C = DC.build(True); M = DC.measured(C)
rows_m = DC.rows
high1 = (rows_m.band.to_numpy()[C.i1] >= 4000); high2 = (rows_m.band.to_numpy()[C.i2] >= 4000)


def within_fn(l0, l1):
    m1 = np.where(C.ice1, pf['q'], np.where(high1, pf['hf'], 1.0)); m2 = np.where(C.ice2, pf['q'], np.where(high2, pf['hf'], 1.0))
    v1 = m1 * pf['a'] * np.exp(-C.z1 / pf['zs']); v2 = m2 * pf['a'] * np.exp(-C.z2 / pf['zs'])
    omr = 1 - (pf['w'] * np.cos(W12 * C.lag) + (1 - pf['w']) * np.exp(-C.lag / pf['tau']))
    dr = 0.5 * (pf['beta'] * C.lag / 12) ** 2
    a, b = l0 + l1 * C.z1.to_numpy(), l0 + l1 * C.z2.to_numpy()
    L2 = 0.5 * (a * a + b * b)
    c = np.sqrt(a * b / L2) * K2.matern(np.abs(C.z1 - C.z2).to_numpy() / np.sqrt(L2), 1.5)
    cd = K2.matern(np.abs(C.z1 - C.z2).to_numpy() / 50.0, 1.5)
    cov = 2 * np.sqrt(v1 * v2) * c * omr + 2 * dr * cd
    return cov, 2 * v1 * omr + 2 * dr, 2 * v2 * omr + 2 * dr


def sse(x):
    l0, l1 = np.exp(x[0]), x[1]
    if l1 < 0:
        return 1e9
    cov, v1, v2 = within_fn(l0, l1)
    t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
    return float(np.sum(M['n'] * (M['corr'] - t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())) ** 2))


best = min((minimize(sse, [np.log(s0), s1], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxfev=4000))
            for s0 in (0.1, 0.5, 2) for s1 in (0.1, 0.5, 1.0)), key=lambda r: r.fun)
l0, l1 = float(np.exp(best.x[0])), float(best.x[1])
print(f'depth length l(z) = {l0:.3f} + {l1:.3f} z  (SSE {best.fun:.1f})')
cov, v1, v2 = within_fn(l0, l1)
t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb)).groupby(['zb', 'dzb'], observed=True)
M['model'] = t.c.sum() / np.sqrt(t.v1.sum() * t.v2.sum())
M['corr_med'] = pd.concat([DC.measured(C, r)['corr'] for r in range(32)], axis=1).median(axis=1)
print(M.round(3).to_string())
FINAL = dict(pf, l0=l0, l1=l1, lb0=l0, lb1=l1, ld=50.0, nu=1.5, nub=1.5, centre=True)
pickle.dump(dict(final=FINAL, corr_table=M), open(W + 'final_kernel.pkl', 'wb'))
print('FINAL', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in FINAL.items()})
