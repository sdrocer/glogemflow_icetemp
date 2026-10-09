import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2lib as L, s2vfit as V
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = ((rows.band >= 4400) & (Y[13] > -4.0) & (y < -8.0) & (rows.depth > 1.5)).to_numpy()
np.save(W + 'spike_rows_month.npy', np.flatnonzero(spike))
P = pd.read_pickle(W + 'pairs.pkl')
keep = ~(P.i.isin(np.flatnonzero(spike)) | P.j.isin(np.flatnonzero(spike))).to_numpy()
P = P[keep].reset_index(drop=True)
nrows = rows.groupby('t2').size()
nr = (rows.t2.map(nrows) - rows.t2.map(rows[~spike].groupby('t2').size()).fillna(0)).to_numpy()  # unused
cnt = rows[~spike].groupby('t2').size()
w = (2.0 / P.t2.map(cnt)).to_numpy()
z, lag, vobs = P.depth.to_numpy(), P.dlag.to_numpy(), P.var_obs.to_numpy()
r = Y[13] - y
d2 = (r[P.i] - r[P.j]) ** 2
C = {
 'A per+noise':               dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0),
 'B per+noise+slow':          dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0, al=0.5, zl=10.0, tl=24.0),
 'C B + deep floor':          dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0, al=0.5, zl=10.0, bl=0.3, tl=24.0),
 'D per+noise+slow+decadal':  dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0, al=0.5, zl=10.0, tl=24.0, ad=1.0, td=300.0),
 'E D, decadal fades in z':   dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0, al=0.5, zl=10.0, tl=24.0, ad=1.0, td=300.0, zd=50.0),
 'F D, quasi-periodic':       dict(ap=0.5, ar=0.3, zs=4.0, tr=1.0, tp=60.0, al=0.5, zl=10.0, tl=24.0, ad=1.0, td=300.0),
}
fits = {}
for name, st in C.items():
    p, cl, k = V.fit(st, z, lag, d2, vobs, w)
    fits[name] = (p, cl, k)
    print(f'{name:26s} CL {cl:9.2f} k {k:2d}  ' + ' '.join(f'{a}={b:.4g}' for a, b in p.items()))
pickle.dump(dict(fits=fits, spike=np.flatnonzero(spike)), open(W + 'vfit_d13.pkl', 'wb'))
# compare binned
P['meas'] = d2 / 2
lagnames = [n for _, _, n in L.LAG_BINS]
for name in ('B per+noise+slow', 'D per+noise+slow+decadal', 'E D, decadal fades in z'):
    P['mod'] = V.gamma_model(fits[name][0], z, lag) + vobs
    t = P.groupby(['zb', 'lagbin']).agg(meas=('meas', 'mean'), model=('mod', 'mean')).unstack('zb')
    print(f'\n== {name} (binned, pair-weighted) =='); print(t.reindex(lagnames).round(3).to_string())
