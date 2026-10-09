"""Baseline fits: T2 and T3r, all common rows; surface, emulator mean, emulator with plug-in E."""
import sys, time, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ctx = C.context()
res = {}
for name, c in ctx.items():
    t0 = time.time()
    r = C.fit_rows(c, np.ones(len(c.y), bool), methods=('surface', 'emulator', 'emulator_E'))
    res[name] = {k: {kk: vv for kk, vv in v.items() if kk not in ('grid', 'J')} for k, v in r.items()}
    for k, s in r.items():
        G = C.predict_at(c, s, 'surface')
        extra = (f" | pass1 MAP {s['pass1_map'][0]:.3f}/{s['pass1_map'][1]:.3f}, E median var {s['E_median_var']:.4f}, "
                 f"top eigenvalue {s['E_top_eig']:.1f}") if k == 'emulator_E' else ''
        print(f'{name} {k:10s}: {C.fmt(s)} | RMSE at MAP (surface) {C.rmse(G, c.y):.3f}{extra}  [{time.time()-t0:.0f} s]')
pickle.dump(res, open(C.OUT / 'baseline.pkl', 'wb'))
