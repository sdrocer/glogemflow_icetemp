"""(e) supplement: the linearised summary (Fisher + unit prior precision, prior-sd units) at the
5 best runs by T2 RMSE (13, 14, 3, 9, 0; f_adv set to 1.0) for T2 and T3r with and without the top
3 m, both methods. The decision rule's 'told apart': f_ins spread <= 0.5 and corr <= 0.8."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
pd.set_option('display.width', 250, 'display.max_columns', 40)
ctx = C.context()
rec = []
for name in ('T2', 'T3r'):
    c = ctx[name]
    for lab, m in (('all rows', np.ones(len(c.y), bool)), ('without top 3 m', (c.rows.depth > 3.0).to_numpy())):
        S = c.Sigma[np.ix_(m, m)]
        for kind in ('surface', 'emulator'):
            sur = c.sur.subset(m) if kind == 'surface' else C.EmuSub(c.emu, np.flatnonzero(m))
            for i in (13, 14, 3, 9, 0):
                th0 = c.d.theta[i].copy(); th0[2] = 1.0
                L = C.linearised(sur, c.y[m], S, th0, None, kind)
                rec.append(dict(target=name, rows=lab, method=kind, run=i, f_ref=th0[0], f_ins=th0[1], corr=L['corr'],
                                spread_f_ref=L['spread_f_ref'], spread_f_ins=L['spread_f_ins'],
                                passes=bool(L['corr'] <= 0.8 and L['spread_f_ins'] <= 0.5)))
D = pd.DataFrame(rec)
print(D.round(3).to_string(index=False))
print(D.groupby(['target', 'rows', 'method']).agg(passes=('passes', 'sum'), corr_max=('corr', 'max'), spread_ins_max=('spread_f_ins', 'max')).round(3).to_string())
