"""Emulator refit noise: the same Y refitted, and Y plus tiny noise (0.01 C), optimum moves."""
import sys
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ctx = C.context()
rng = np.random.default_rng(0)
for name in ('T2', 'T3r'):
    c = ctx[name]
    out = []
    for rep in range(6):
        Y = c.Y.copy()
        if rep >= 2:
            Y = Y + rng.normal(0, 0.01, Y.shape)
        e = C.Emu(c.d.theta[c.d.keep], Y[c.d.keep], 0.999, c.cap)
        c2 = C.SimpleNamespace(**vars(c)); c2.emu = e
        r = C.fit_rows(c2, np.ones(len(c.y), bool), methods=('emulator',))['emulator']
        out.append((r['mean_f_ref'], r['mean_f_ins'], r['corr'], r['spread_f_ins']))
        print(f'{name} rep {rep} ({"same Y" if rep < 2 else "Y + N(0, 0.01 C)"}): p {e.p} mean f_ref {r["mean_f_ref"]:.3f} f_ins {r["mean_f_ins"]:.3f} corr {r["corr"]:+.2f} spread f_ins {r["spread_f_ins"]:.2f}')
    o = np.array(out)
    print(f'{name}: spread over refits: f_ref sd {o[:,0].std():.4f} ({o[:,0].std()/C.PSD[0]:.3f} prior sd), f_ins sd {o[:,1].std():.4f} ({o[:,1].std()/C.PSD[1]:.3f} prior sd)')
