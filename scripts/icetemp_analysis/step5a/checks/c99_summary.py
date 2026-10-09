"""Collect the Stage 3c key numbers from the check outputs into c99_summary.json."""
import sys, json, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
O = C.OUT
keep = ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref', 'spread_f_ins',
        'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd', 'mahal_base', 'chi2_n']
out = {'baseline': {t: {k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk in
                            ('map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'sd_f_ref', 'sd_f_ins', 'corr', 'spread_f_ref',
                             'spread_f_ins', 'lin_corr', 'lin_spread_f_ref', 'lin_spread_f_ins', 'ratio_map', 'E_top_eig')}
                        for k, v in r.items()} for t, r in pickle.load(open(O / 'baseline.pkl', 'rb')).items()}}
for key, f in (('a_firn_ice', 'c02_firn_ice.pkl'), ('b_hef', 'c03_hef.pkl'), ('c_melting', 'c05_melting.pkl'),
               ('d_snow', 'c04c_snow_fit.pkl'), ('e_f_top3m_cg', 'c06_top3m_cg.pkl'), ('g_old_new', 'c07_old_new.pkl')):
    p = pickle.load(open(O / f, 'rb'))
    T = p['table']
    cols = [c for c in T.columns if c in keep or c.startswith('rmse')]
    out[key] = json.loads(T[cols].round(4).to_json(orient='records'))
    if key == 'a_firn_ice':
        out['a_z'] = json.loads(p['z'].round(4).to_json(orient='records'))
        out['a_groups'] = json.loads(p['groups'].round(4).to_json(orient='records'))
    if key == 'g_old_new':
        out['g_residuals'] = json.loads(p['residuals'].round(4).to_json(orient='records'))
json.dump(out, open(O / 'c99_summary.json', 'w'), indent=1)
print('written', O / 'c99_summary.json', 'sections', list(out))
