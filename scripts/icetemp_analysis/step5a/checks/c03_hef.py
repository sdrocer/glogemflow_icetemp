"""(b) The fit with and without the Hintereisferner tongue (ice, bands < 3000 m: 2555, 2615, 2805,
2955 m), T2 and T3r, f_adv = 1, two methods. Writes c03_hef.log / .pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
pd.set_option('display.width', 280, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
HEF = 'Hintereisferner tongue, ice (< 3000 m)'
c2, c3 = ctx['T2'], ctx['T3r']
h2, h3 = (c2.rows.grp == HEF).to_numpy(), (c3.rows.grp == HEF).to_numpy()
print(f'HEF tongue rows: T2 {h2.sum()} of {len(h2)}, T3r {h3.sum()} of {len(h3)}')
masks2 = {'T2 all': np.ones(len(c2.y), bool), 'T2 without HEF tongue': ~h2, 'HEF tongue': h2}
variants = [('all rows', 'T2', np.ones(len(c2.y), bool), None), ('without HEF tongue', 'T2', ~h2, None),
            ('all rows', 'T3r', np.ones(len(c3.y), bool), None), ('without HEF tongue', 'T3r', ~h3, None)]
recs, fits = C.run_variants(ctx, variants, base, masks2)
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref',
                   'spread_f_ins', 'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd', 'd_f_ref_base_sd',
                   'd_f_ins_base_sd', 'mahal_base', 'chi2_n', 'rmse T2 all', 'rmse T2 without HEF tongue', 'rmse HEF tongue'])
print(T.round(3).to_string(index=False))
# HEF tongue residual at both optima (bias), and at design 13
for kind in ('surface', 'emulator'):
    for lab in ('all rows', 'without HEF tongue'):
        s = fits[(lab, 'T2')][kind]
        th = np.array([s['map_f_ref'], s['map_f_ins'], 1.0])
        G = c2.sur.predict(th)[0] if kind == 'surface' else c2.emu.predict(th)[0][0]
        print(f'{kind:8s} T2 fit {lab:20s}: HEF tongue bias {np.mean(G[h2] - c2.y[h2]):+.3f} C, rest bias {np.mean(G[~h2] - c2.y[~h2]):+.3f} C')
print(f'design 13 (actual run): HEF tongue bias {np.mean(c2.Y[13, h2] - c2.y[h2]):+.3f} C, RMSE {C.rmse(c2.Y[13, h2], c2.y[h2]):.3f}')
pickle.dump(dict(table=T, fits=fits), open(C.OUT / 'c03_hef.pkl', 'wb'))
