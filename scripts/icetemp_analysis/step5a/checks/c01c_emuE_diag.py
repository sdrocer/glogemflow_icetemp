"""Why the emulator with its own uncertainty E is used for reference only: E evaluated at three
points of the T2 grid (the emulator optimum, the surface optimum, design 13), with the raw and
the leave-one-run-out inflated GP variance; the T2 posterior with each E plugged in."""
import sys
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ctx = C.context(('T2',))
c = ctx['T2']
m = np.ones(len(c.y), bool)
sur = C.EmuSub(c.emu, np.flatnonzero(m))
trust = C.trust_mask_fn(c.d.theta[c.d.keep])
for infl in (1.0, 1.86):
    c.emu.infl = infl
    for th in ([0.311, 0.659, 1.0], [0.333, 0.762, 1.0], [0.524, 1.089, 1.0]):
        E = C.emu_cov(sur, th)
        ev = np.linalg.eigvalsh(E)[::-1][:4]
        s = C.fit(sur, c.y, c.Sigma + E, None, 'emulator', trust)
        print(f'infl {infl} E at {th}: top eigenvalues {np.round(ev, 2)} median diag {np.median(np.diag(E)):.4f} -> '
              f'MAP {s["map_f_ref"]:.3f}/{s["map_f_ins"]:.3f} mean {s["mean_f_ref"]:.3f}/{s["mean_f_ins"]:.3f} corr {s["corr"]:+.2f} spread {s["spread_f_ref"]:.2f}/{s["spread_f_ins"]:.2f}')
