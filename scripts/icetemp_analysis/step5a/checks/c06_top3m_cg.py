"""(e) T3 (= T3r, Stage 2) without the top 3 m (point depth <= 3.0 m left out; logger depths drift
with burial or melt), with T2 for comparison; (f) Colle Gnifetti deeper than 40 m (Monte Rosa
entities at band >= 4400 m, point depth > 40 m: the accepted pre-1942 limitation) in vs out.
f_adv = 1, two methods. Writes c06_top3m_cg.log / .pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
ta = C.ta
pd.set_option('display.width', 280, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
c2, c3 = ctx['T2'], ctx['T3r']
CG = 'Colle Gnifetti (>= 4400 m)'
top2, top3 = (c2.rows.depth <= 3.0).to_numpy(), (c3.rows.depth <= 3.0).to_numpy()
cgd2 = ((c2.rows.grp == CG) & (c2.rows.depth > 40.0)).to_numpy()
cgd3 = ((c3.rows.grp == CG) & (c3.rows.depth > 40.0)).to_numpy()
for c, top, cgd in ((c2, top2, cgd2), (c3, top3, cgd3)):
    print(f'{c.name}: top 3 m rows {top.sum()} of {len(top)}; by group {c.rows[top].grp.value_counts().to_dict()}')
    print(f'{c.name}: Colle Gnifetti deeper than 40 m: {cgd.sum()} rows, entities {sorted(c.rows[cgd].entity.unique())}, '
          f'depth {c.rows.depth[cgd].min():.1f}-{c.rows.depth[cgd].max():.1f} m; design 13 bias {np.mean(c.Y[13, cgd] - c.y[cgd]):+.2f} C')
# logger entities (many months at one depth): rows per point
for c in (c3,):
    npt = c.rows.groupby('t2').size()
    lg = c.rows[c.rows.t2.map(npt) >= 6]
    print('T3r points with >= 6 rows (logger-like):', lg.t2.nunique(), 'points in', sorted(lg.entity.unique()))
    print('   of them in the top 3 m:', lg[lg.depth <= 3].t2.nunique(), 'points,', int((lg.depth <= 3).sum()), 'rows')
masks2 = {'T2 all': np.ones(len(c2.y), bool), 'T2 top 3 m': top2, 'T2 below 3 m': ~top2, 'CG deep': cgd2, 'T2 without CG deep': ~cgd2}
variants = [('all rows', 'T2', np.ones(len(c2.y), bool), None), ('without top 3 m', 'T2', ~top2, None),
            ('all rows', 'T3r', np.ones(len(c3.y), bool), None), ('without top 3 m', 'T3r', ~top3, None),
            ('without CG > 40 m', 'T2', ~cgd2, None), ('without CG > 40 m', 'T3r', ~cgd3, None),
            ('without top 3 m and CG > 40 m', 'T2', ~top2 & ~cgd2, None), ('without top 3 m and CG > 40 m', 'T3r', ~top3 & ~cgd3, None)]
recs, fits = C.run_variants(ctx, variants, base, masks2, methods=('surface', 'emulator', 'emulator_E'))
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref',
                   'spread_f_ins', 'lin_corr', 'lin_spread_f_ref', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd',
                   'd_f_ref_base_sd', 'd_f_ins_base_sd', 'mahal_base', 'chi2_n', 'rmse T2 all', 'rmse T2 top 3 m', 'rmse CG deep',
                   'rmse T2 without CG deep'])
print(T.round(3).to_string(index=False))
pickle.dump(dict(table=T, fits=fits), open(C.OUT / 'c06_top3m_cg.pkl', 'wb'))
