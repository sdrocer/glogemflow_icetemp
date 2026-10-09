"""(a) Do firn points and ice points want the same f_ins? T2 and T3r, f_adv = 1, two methods.
Splits: the borehole's firn label (data.py: firn / ice / unknown) and the model band's state at
the readings' dates (firnice_firnline: ice band = the snow-on-ice branch where f_ins scales the
snow conductivity; firn band = f_ins scales the whole firn column). Also every borehole group's
own f_ins with f_ref held at the joint optimum. Writes c02_firn_ice.log / .pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
ta = C.ta
pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
SR = pickle.load(open(C.OUT / 'snow_rows.pkl', 'rb'))['S']           # per T3 row: ice_m (model band ice in that month)
obs = ctx['T2'].d.obs
# model ice-band share per T2 point (count-weighted over its months) and per T3r row (mean over its months)
ice_t3 = dict(zip(zip(SR.t2, SR.ym), SR.ice_m.astype(float)))
w = SR.assign(wi=SR['count'] * SR.ice_m).groupby('t2')[['wi', 'count']].sum()
ice_t2 = (w.wi / w['count']).to_dict()
c2, c3 = ctx['T2'], ctx['T3r']
c2.rows['model_ice'] = c2.rows.t2.map(ice_t2).to_numpy()
c3.rows['model_ice'] = [np.mean([ice_t3[(t, int(mm))] for mm in ms]) for t, ms in zip(c3.rows.t2, c3.rows.months)]
for c in (c2, c3):
    print(f'{c.name}: model ice-band share per row: 0 -> {int((c.rows.model_ice == 0).sum())}, 1 -> {int((c.rows.model_ice == 1).sum())}, '
          f'between -> {int(((c.rows.model_ice > 0) & (c.rows.model_ice < 1)).sum())}')
    print(pd.crosstab(c.rows.firn_label, c.rows.model_ice >= 0.5, margins=True).rename(columns={False: 'model firn band', True: 'model ice band'}).to_string())
    print(pd.crosstab(c.rows.grp, [c.rows.firn_label, c.rows.model_ice >= 0.5]).to_string())

base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
recs, out = [], {}
sets = lambda c: [('all rows', np.ones(len(c.y), bool)),
                  ('label firn', (c.rows.firn_label == 'firn').to_numpy()),
                  ('label ice', (c.rows.firn_label == 'ice').to_numpy()),
                  ('label unknown', (c.rows.firn_label == 'unknown').to_numpy()),
                  ('label firn + unknown', (c.rows.firn_label != 'ice').to_numpy()),
                  ('model firn band', (c.rows.model_ice < 0.5).to_numpy()),
                  ('model ice band', (c.rows.model_ice >= 0.5).to_numpy()),
                  ('model ice band, top 15 m', ((c.rows.model_ice >= 0.5) & (c.rows.depth <= 15)).to_numpy()),
                  ('model ice band, without Hintereisferner tongue', ((c.rows.model_ice >= 0.5) & (c.rows.grp != 'Hintereisferner tongue, ice (< 3000 m)')).to_numpy())]
for c in (c2, c3):
    for lab, m in sets(c):
        r = C.fit_rows(c, m)
        out[(c.name, lab)] = C.strip(r)
        recs += C.row(lab, c.name, r, base[c.name])
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'sd_f_ref', 'mean_f_ins', 'sd_f_ins',
                   'corr', 'spread_f_ref', 'spread_f_ins', 'd_f_ins_prior_sd', 'chi2_n'])
print('\n== (a) fits on firn and ice rows (f_adv = 1; grid posterior; spreads = posterior sd / prior sd) ==')
print(T.round(3).to_string(index=False))
# z-scores firn vs ice on f_ins (independent row sets)
print('\n== f_ins wanted by firn rows vs ice rows: difference / sqrt(sd_a^2 + sd_b^2) ==')
zrec = []
for tname in ('T2', 'T3r'):
    for a, b in (('label firn', 'label ice'), ('label firn + unknown', 'label ice'), ('model firn band', 'model ice band'),
                 ('model firn band', 'model ice band, top 15 m'), ('model firn band', 'model ice band, without Hintereisferner tongue')):
        for kind in ('surface', 'emulator'):
            sa, sb = out[(tname, a)][kind], out[(tname, b)][kind]
            dz = (sb['mean_f_ins'] - sa['mean_f_ins']) / np.hypot(sa['sd_f_ins'], sb['sd_f_ins'])
            zrec.append(dict(target=tname, a=a, b=b, method=kind, f_ins_a=sa['mean_f_ins'], sd_a=sa['sd_f_ins'], f_ins_b=sb['mean_f_ins'],
                             sd_b=sb['sd_f_ins'], diff=sb['mean_f_ins'] - sa['mean_f_ins'],
                             diff_prior_sd=(sb['mean_f_ins'] - sa['mean_f_ins']) / C.PSD[1], z=dz))
Z = pd.DataFrame(zrec)
print(Z.round(3).to_string(index=False))
# each borehole group's f_ins with f_ref at the joint optimum (T2 surface MAP)
print('\n== each group alone: f_ins posterior with f_ref fixed at the joint T2 / T3r optimum (prior included); data-only optimum ==')
grec = []
for c in (c2, c3):
    for kind in ('surface', 'emulator'):
        fr = base[c.name][kind]['map_f_ref']
        for g in sorted(c.rows.grp.unique()):
            for lab, m in ((g, (c.rows.grp == g).to_numpy()),):
                s = C.fit_fins_given_fref(c, m, fr, kind=kind)
                G = C.predict_at(c, base[c.name][kind], 'surface')[m]
                grec.append(dict(target=c.name, method=kind, group=g, n=s['n'], f_ref_fixed=fr, f_ins_mean=s['mean'], f_ins_sd=s['sd'],
                                 f_ins_map=s['map'], f_ins_data_opt=s['data_opt'], at_bound=s['at_bound'],
                                 dchi2_over_range=s['dchi2_range'], bias_at_joint=float(np.mean(G - c.y[m])),
                                 label=c.rows.firn_label[m].mode().iat[0], model_ice=float(c.rows.model_ice[m].mean())))
Gt = pd.DataFrame(grec)
print(Gt.round(3).to_string(index=False))
pickle.dump(dict(fits=out, table=T, z=Z, groups=Gt), open(C.OUT / 'c02_firn_ice.pkl', 'wb'))
