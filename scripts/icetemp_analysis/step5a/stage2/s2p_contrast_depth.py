"""Contrast part of the T3r misfit by depth band and system: per-row r_i (Pi r)_i with
Pi = Sigma^-1 - Sigma^-1 Z (Z' Sigma^-1 Z)^-1 Z' Sigma^-1 (Z: pooled points); dims = rows - points."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve
import target_analysis as ta, s2kern2 as K2
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
FINAL = pickle.load(open(W + 'final_kernel.pkl', 'rb'))['final']
t = ta.target('T3r', obs, mod)
rows, r = t.rows, t.Y[13] - t.y
Sig = ta.error_cov(rows) + K2.seasonal_matrix(rows, FINAL, obs=obs) + np.diag(K2.date_spread(rows, obs, mod['Y3'][13]))
cc = np.zeros(len(r))
for ent in np.unique(rows.e):
    idx = np.flatnonzero(rows.e.to_numpy() == ent)
    c = cho_factor(Sig[np.ix_(idx, idx)], lower=True)
    pts, inv = np.unique(rows.t2.to_numpy()[idx], return_inverse=True)
    Z = np.zeros((len(idx), len(pts))); Z[np.arange(len(idx)), inv] = 1
    SiZ = cho_solve(c, Z); Sir = cho_solve(c, r[idx])
    Pr = Sir - SiZ @ np.linalg.solve(Z.T @ SiZ, Z.T @ Sir)
    cc[idx] = r[idx] * Pr
rows = rows.assign(cc=cc, zb=ta.depth_band(rows.depth))
npt = rows.groupby('t2').size()
rows['dim'] = 1 - 1 / rows.t2.map(npt)          # each row carries (n-1)/n of a contrast dimension
g = rows.groupby(['zb']).agg(chi2=('cc', 'sum'), dims=('dim', 'sum'))
g['ratio'] = g.chi2 / g.dims
print(g.round(3).to_string())
g = rows.groupby(['system_short', 'zb']).agg(chi2=('cc', 'sum'), dims=('dim', 'sum'))
g['ratio'] = g.chi2 / g.dims
print(g[g.dims > 3].round(3).to_string())
