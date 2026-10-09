"""(g) supplement: T3r misfit ratio of old and new rows split into level and month-to-month parts
(target_analysis.misfit_ratio) at the T3r surface optimum, per entity with old rows."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ta = C.ta
ctx = C.context(('T3r',))
c = ctx['T3r']
s = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))['T3r']['surface']
G = c.sur.predict(np.array([s['map_f_ref'], s['map_f_ins'], 1.0]))[0]
r = G - c.y
old = (c.rows.ym_last < 197901).to_numpy(); new = (c.rows.ym_first >= 197901).to_numpy()
HEF = 'Hintereisferner tongue, ice (< 3000 m)'
for lab, m in (('old', old), ('new', new), ('old, Hintereisferner tongue', old & (c.rows.grp == HEF).to_numpy()),
               ('old, rest', old & (c.rows.grp != HEF).to_numpy())):
    mr = ta.misfit_ratio(r[m], c.Sigma[np.ix_(m, m)], c.rows[m].reset_index(drop=True))
    print(f'{lab:30s}: rows {mr["n"]}, points {mr["points"]}, ratio {mr["ratio"]:.2f}, level {mr["level"]:.2f}, month-to-month {mr["month_to_month"]:.2f}')
for e in sorted(c.rows[old].entity.unique()):
    m = old & (c.rows.entity == e).to_numpy()
    mr = ta.misfit_ratio(r[m], c.Sigma[np.ix_(m, m)], c.rows[m].reset_index(drop=True))
    print(f'   {e:28s}: rows {mr["n"]:3d}, points {mr["points"]:3d}, ratio {mr["ratio"]:.2f}, level {mr["level"]:.2f}, month-to-month {mr["month_to_month"]:.2f}')
