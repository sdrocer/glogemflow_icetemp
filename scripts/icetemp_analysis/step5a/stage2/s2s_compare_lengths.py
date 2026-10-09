"""T3r / T3m / T3 misfit with linear, knot and monotone-knot depth lengths; contrast part by depth."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
import target_analysis as ta, s2kern2 as K2
import s2r_length_knots as KN
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
d = pickle.load(open(W + 'final_kernel.pkl', 'rb'))
# monotone: one length from 20 m down
best = min((minimize(lambda x: KN.sse(np.log([np.exp(x[0]), np.exp(x[1]), np.exp(x[2]), np.exp(x[2])])), np.log(s0),
                     method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxfev=6000))
            for s0 in ([1.2, 4.5, 60], [1.2, 4.5, 300], [1.0, 3.0, 20])), key=lambda r: r.fun)
l3 = np.exp(best.x)
lk_m = [l3[0], l3[1], l3[2], l3[2]]
print('monotone knot lengths (2.5, 9, >=20 m):', np.round(l3, 2), f'SSE {best.fun:.1f}')
d['final_mono'] = dict(d['final'], lk=lk_m, lkb=lk_m)
pickle.dump(d, open(W + 'final_kernel.pkl', 'wb'))


def contrast_by_depth(rows, r, Sig):
    cc = np.zeros(len(r))
    for ent in np.unique(rows.e):
        idx = np.flatnonzero(rows.e.to_numpy() == ent)
        c = cho_factor(Sig[np.ix_(idx, idx)], lower=True)
        pts, inv = np.unique(rows.t2.to_numpy()[idx], return_inverse=True)
        Z = np.zeros((len(idx), len(pts))); Z[np.arange(len(idx)), inv] = 1
        SiZ = cho_solve(c, Z); Sir = cho_solve(c, r[idx])
        cc[idx] = r[idx] * (Sir - SiZ @ np.linalg.solve(Z.T @ SiZ, Z.T @ Sir))
    npt = rows.groupby('t2').size()
    t = pd.DataFrame(dict(cc=cc, dim=1 - 1 / rows.t2.map(npt).to_numpy(), zb=ta.depth_band(rows.depth)))
    g = t.groupby('zb').agg(chi2=('cc', 'sum'), dims=('dim', 'sum'))
    return (g.chi2 / g.dims).round(3).to_dict()


out = []
for lab in ('final', 'final_knots', 'final_mono'):
    p = d[lab]
    for name in ('T3r', 'T3m', 'T3'):
        t = ta.target('T3' if name == 'T3m' else name, obs, mod)
        rr, Yt, yt = t.rows, t.Y, t.y
        if name == 'T3m':
            m = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[m].reset_index(drop=True), Yt[:, m], yt[m]
        Sig = ta.error_cov(rr) + K2.seasonal_matrix(rr, p, obs=obs)
        if name == 'T3r':
            Sig = Sig + np.diag(K2.date_spread(rr, obs, mod['Y3'][13]))
        r = Yt[13] - yt
        mm = K2.misfit(r, Sig, rr)
        rec = dict(lengths=lab, target=name, ratio=mm['ratio'], level=mm['level'], contrast=mm['contrast'])
        if name == 'T3r':
            rec.update({f'contrast {k}': v for k, v in contrast_by_depth(rr, r, Sig).items()})
        out.append(rec)
print(pd.DataFrame(out).round(3).to_string(index=False))
