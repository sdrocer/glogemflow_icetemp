"""(d) The depth reference under snow: how many shallow (<= 3 m) T2 points / T3r rows were read
while the model's ice band carried snow, and how the fit changes if they are left out or read
z + h deep (h = modelled snow depth; c04b_snow_shift.py): only the shallow rows, or every depth.
T2 and T3r, f_adv = 1, two methods. Writes c04c_snow_fit.log / .pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
ta = C.ta
pd.set_option('display.width', 280, 'display.max_columns', 40, 'display.max_rows', 200)
SR = pickle.load(open(C.OUT / 'snow_rows.pkl', 'rb'))['S']
SH = pickle.load(open(C.OUT / 'snow_shift.pkl', 'rb'))
d = C.load()
mod, obs = d.mod, d.obs
c2m, c3m = mod['common2'], mod['common3r']
snow = dict(zip(zip(SR.t2, SR.ym), SR.snow_ice))
hm = dict(zip(zip(SR.t2, SR.ym), SR.h_m))
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
ctx0 = C.context()
r2, r3 = ctx0['T2'].rows, ctx0['T3r'].rows
# flags
r3_snow = np.array([any(snow[(t, int(m))] for m in ms) for t, ms in zip(r3.t2, r3.months)])
sh3 = (r3.depth <= 3.0).to_numpy()
pt_snow = SR[SR.snow_ice & (SR.depth <= 3.0)].t2.unique()
sh2 = (r2.depth <= 3.0).to_numpy()
p2_snow = r2.t2.isin(pt_snow).to_numpy()
w_all = SR.groupby('t2')['count'].sum(); w_sn = SR[SR.snow_ice].groupby('t2')['count'].sum()
print(f'T2 points <= 3 m: {sh2.sum()}; read (in part) while the model ice band carried snow: {p2_snow.sum()} '
      f'({w_sn.reindex(r2.t2[p2_snow]).fillna(0).sum():.0f} of their {w_all.reindex(r2.t2[p2_snow]).sum():.0f} readings)')
print(r2[p2_snow].groupby('entity').size().to_dict())
print(f'T3r rows <= 3 m: {sh3.sum()}; under snow on an ice band: {(sh3 & r3_snow).sum()}; at any depth: {r3_snow.sum()}')
print(r3[sh3 & r3_snow].groupby('entity').size().to_dict())
hh = np.array([np.mean([hm[(t, int(m))] for m in ms]) for t, ms in zip(r3.t2, r3.months)])
print(f'modelled snow depth at those shallow rows: median {np.median(hh[sh3 & r3_snow]):.2f} m, range {hh[sh3 & r3_snow].min():.2f}-{hh[sh3 & r3_snow].max():.2f} m')
# change of the model values (design 13) at those rows
for v in ('shallow', 'all'):
    dY3r = SH[v]['Y3r'][:, c3m] - mod['Y3r'][:, c3m]
    dY2 = SH[v]['Y2'][:, c2m] - mod['Y2'][:, c2m]
    ch3 = np.abs(dY3r).max(0) > 1e-9
    ch2 = np.abs(dY2).max(0) > 1e-9
    print(f'{v}: T3r rows changed {ch3.sum()}, design 13 change mean {dY3r[13, ch3].mean():+.3f} C, RMS {np.sqrt((dY3r[13, ch3] ** 2).mean()):.3f}; '
          f'T2 points changed {ch2.sum()}, mean {dY2[13, ch2].mean():+.3f}, RMS {np.sqrt((dY2[13, ch2] ** 2).mean()):.3f}; '
          f'misfit at those T3r rows, design 13: before {np.mean(mod["Y3r"][13, c3m][ch3] - r3.obs[ch3]):+.3f} after {np.mean(SH[v]["Y3r"][13, c3m][ch3] - r3.obs[ch3]):+.3f}')
    assert np.isfinite(SH[v]['Y2'][:, c2m]).all() and np.isfinite(SH[v]['Y3r'][:, c3m]).all()
masks2 = {'T2 all': np.ones(len(r2), bool), 'T2 shallow under snow': p2_snow}
recs, fits = C.run_variants(ctx0, [('as now', 'T2', np.ones(len(r2), bool), None),
                                   ('shallow under snow left out', 'T2', ~p2_snow, None),
                                   ('as now', 'T3r', np.ones(len(r3), bool), None),
                                   ('shallow under snow left out', 'T3r', ~(sh3 & r3_snow), None)], base, masks2)
for v in ('shallow', 'all'):
    ctxv = C.context(Y_override={'T2': SH[v]['Y2'][:, c2m], 'T3r': SH[v]['Y3r'][:, c3m]})
    rr, ff = C.run_variants(ctxv, [(f'depth + snow, {v} rows', 'T2', np.ones(len(r2), bool), None),
                                   (f'depth + snow, {v} rows', 'T3r', np.ones(len(r3), bool), None)], base, masks2)
    recs += rr; fits.update(ff)
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref',
                   'spread_f_ins', 'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd', 'd_f_ref_base_sd',
                   'd_f_ins_base_sd', 'mahal_base', 'chi2_n', 'rmse T2 all', 'rmse T2 shallow under snow'])
print(T.round(3).to_string(index=False))
pickle.dump(dict(table=T, fits=fits, p2_snow=p2_snow, r3_snow_shallow=sh3 & r3_snow), open(C.OUT / 'c04c_snow_fit.pkl', 'wb'))
