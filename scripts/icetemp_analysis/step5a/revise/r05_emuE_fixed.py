"""'Told apart' with the emulator's uncertainty E held fixed (revise stage): Stage 3c's emulator
(checks/c3lib.py) with E computed once at a point and kept, for T2 and T3r, at three points each
(Stage 3c's emulator optimum without E, the surface optimum, design 13) and three inflations of E
(1.0; T3r's own leave-one-run-out factor 0.708; T2's 1.86). Exact grid fits as in Stage 3c,
campaign prior. Extends c01c_emuE_diag (T2 only) to T3r. Read-only use of c3lib.
Output: emuE_fixed.json."""
import json
import sys
import time

sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np  # noqa: E402
import c3lib as C  # noqa: E402

t0 = time.time()
ctx = C.context(('T3r', 'T2'))
POINTS = {'T3r': {"Stage 3c's emulator optimum without E": [0.239, 0.609, 1.0], 'surface optimum': [0.273, 0.714, 1.0],
                  'design 13': [0.524, 1.089, 1.0]},
          'T2': {"Stage 3c's emulator optimum without E": [0.311, 0.659, 1.0], 'surface optimum': [0.333, 0.762, 1.0],
                 'design 13': [0.524, 1.089, 1.0]}}
out = []
for name, pts in POINTS.items():
    c = ctx[name]
    sur = C.EmuSub(c.emu, np.flatnonzero(np.ones(len(c.y), bool)))
    trust = C.trust_mask_fn(c.d.theta[c.d.keep])
    for infl in (1.0, 0.708, 1.856):
        c.emu.infl = infl
        for where, th in pts.items():
            E = C.emu_cov(sur, th)
            ev = np.linalg.eigvalsh(E)[::-1]
            s = C.fit(sur, c.y, c.Sigma + E, None, 'emulator', trust)
            ok = abs(s['corr']) <= 0.8 and s['spread_f_ins'] <= 0.5
            out.append(dict(target=name, E_at=where, theta=th[:2], inflation=infl, E_top=float(ev[0]),
                            map_f_ref=float(s['map_f_ref']), map_f_ins=float(s['map_f_ins']), corr=float(s['corr']),
                            spread_f_ref=float(s['spread_f_ref']), spread_f_ins=float(s['spread_f_ins']), told_apart=bool(ok)))
            print(name, where, infl, round(ev[0], 2), round(s['corr'], 3), round(s['spread_f_ins'], 3), ok,
                  f'{time.time() - t0:.0f} s', flush=True)
with open('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/revise/emuE_fixed.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('done', f'{time.time() - t0:.0f} s')
