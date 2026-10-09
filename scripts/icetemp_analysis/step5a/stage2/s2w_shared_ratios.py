"""Misfit ratios with D shared by the months of one borehole at one point (+ n0), S, and the term
refitted for it (median over the 29 runs); compared with the diagonal-D proposal."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2lib as L
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
T = pd.read_pickle(W + 'shared_fit_runs.pkl')
keep = mod['design'].keep.to_numpy()
med = T[keep].median()
PS = dict(a=round(med.a, 2), zs=round(med.zs, 2), ice=round(med.q, 2), high=round(med.hf, 2), w=round(med.w, 2),
          tau=round(med.tau, 1), u0=round(med.u0, 2), beta=round(med.beta, 3))
N0 = 1e-5
print('shared-D term (median29, rounded):', PS)


def obs_cov_shared(rows, n0=N0):
    """sigma_i sigma_j x (borehole overlap) for rows of one point; n0 on the diagonal."""
    sig = rows['sigma'].to_numpy(float); t2 = rows['t2'].to_numpy(int)
    labs = [str(b).split('/') for b in rows['bh_label']]
    D = np.diag(sig ** 2 + n0)
    for pt in np.unique(t2):
        j = np.flatnonzero(t2 == pt)
        if len(j) < 2:
            continue
        names = sorted({x for i in j for x in labs[i]})
        B = np.zeros((len(j), len(names)))
        for a, i in enumerate(j):
            for x in labs[i]:
                B[a, names.index(x)] += 1.0 / len(labs[i])
        B /= np.linalg.norm(B, axis=1, keepdims=True)
        D[np.ix_(j, j)] = np.outer(sig[j], sig[j]) * (B @ B.T) + n0 * np.eye(len(j))
    return D


y2 = ta.target('T2', obs, mod); rm2 = np.sqrt(np.mean((y2.Y - y2.y) ** 2, axis=1)); order = np.argsort(rm2)
runs = [int(x) for x in order[:5]] + [int(order[len(order) // 2])]
out = []
for name in ('T2', 'T3', 'T3r', 'T3 month rows'):
    t = ta.target('T3' if name == 'T3 month rows' else name, obs, mod)
    rr, Yt, yt = t.rows, t.Y, t.y
    if name == 'T3 month rows':
        k = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[k].reset_index(drop=True), Yt[:, k], yt[k]
    S = ta.s_model_cov(rr)
    Dsh = obs_cov_shared(rr)
    Kp = ta.seasonal_cov(rr)                 # diagonal-D proposal (module SEASONAL)
    Ks = ta.seasonal_cov(rr, PS)             # shared-D refit
    for run in runs:
        ds = np.diag(ta.date_spread_var(rr, obs, mod['Y3'][run])) if name == 'T3r' else 0.0
        r = Yt[run] - yt
        a = ta.misfit_ratio(r, np.diag(ta.obs_var(rr)) + S + Kp + ds, rr)
        b = ta.misfit_ratio(r, Dsh + S + Ks + ds, rr)
        out.append(dict(target=name, run=run, diagD_ratio=a['ratio'], diagD_m2m=a['month_to_month'],
                        sharedD_ratio=b['ratio'], sharedD_level=b['level'], sharedD_m2m=b['month_to_month']))
O = pd.DataFrame(out)
print(O.round(3).to_string(index=False))
# T3r by depth and system, design 13, shared D
t = ta.target('T3r', obs, mod); rr = t.rows; r = t.Y[13] - t.y
C = obs_cov_shared(rr) + ta.s_model_cov(rr) + ta.seasonal_cov(rr, PS) + np.diag(ta.date_spread_var(rr, obs, mod['Y3'][13]))
by = []
for lab, mask in ([('depth ' + b, ta.depth_band(rr.depth) == b) for b in L.BANDS]
                  + [('system ' + s, (rr.system_short == s).to_numpy()) for s in sorted(rr.system_short.unique())]):
    mm = ta.misfit_ratio(r[mask], C[np.ix_(mask, mask)], rr[mask].reset_index(drop=True))
    by.append(dict(group=lab, rows=mm['n'], points=mm['points'], ratio=mm['ratio'], level=mm['level'], m2m=mm['month_to_month']))
print(pd.DataFrame(by).round(3).to_string(index=False))
pickle.dump(dict(PS=PS, N0=N0, ratios=O, by=pd.DataFrame(by)), open(W + 'shared_ratios.pkl', 'wb'))
