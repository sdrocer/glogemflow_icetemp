"""Build and check the two surrogates (surface, emulator) for T2 and T3r on the 29 runs with
f_ref >= 0.1; leave-one-run-out accuracy; emulator variance inflation (flat, as the pipeline's
variance_inflation_). Saves surrogates.pkl."""
import sys, time, pickle, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
d = C.load()
keep, theta = d.keep, d.theta
out = {}
for name, t in (('T2', d.t2), ('T3r', d.t3)):
    Y = t.Y[keep]
    cap = C.pmp(t.rows.depth.to_numpy())
    t0 = time.time()
    sur = C.Surface(theta[keep], Y, cap)
    emu = C.Emu(theta[keep], Y, 0.999, cap)
    # leave one run out for the emulator
    idx = np.flatnonzero(keep)
    loo_e = np.zeros_like(Y); z2 = []
    for j in range(len(idx)):
        m = np.arange(len(idx)) != j
        e = C.Emu(theta[keep][m], Y[m], 0.999, cap)
        mean, cvar = e.predict(theta[keep][j])
        loo_e[j] = mean[0] - Y[j]
        var = (e.basis ** 2) @ cvar[0] + e.trunc
        z2.append(loo_e[j] ** 2 / var)
    z2 = np.array(z2)
    infl = float(np.mean(z2))
    emu.infl = infl
    loo_s = sur.loo
    # cap effect on surface loo: predict held out with cap
    print(f'{name}: rows {Y.shape[1]}, emulator p = {emu.p}, trunc var {emu.trunc:.2e}; '
          f'LOO RMSE over all rows/runs: surface {np.sqrt((loo_s**2).mean()):.3f} (uncapped), emulator {np.sqrt((loo_e**2).mean()):.3f}; '
          f'per-row LOO RMSE median surface {np.median(sur.loo_rmse):.3f} p90 {np.percentile(sur.loo_rmse,90):.3f} max {sur.loo_rmse.max():.3f}; '
          f'emulator median {np.median(np.sqrt((loo_e**2).mean(0))):.3f} p90 {np.percentile(np.sqrt((loo_e**2).mean(0)),90):.3f}; '
          f'emulator z^2 mean (inflation) {infl:.2f}, median {np.median(z2):.2f} [{time.time()-t0:.0f} s]')
    # surface vs actual at design 13
    y = t.y
    print(f'   design 13: actual RMSE {C.rmse(t.Y[13], y):.4f}; surface at its theta {C.rmse(sur.predict(theta[13])[0], y):.4f}; '
          f'emulator {C.rmse(emu.predict(theta[13])[0][0], y):.4f}')
    # by group
    rows = t.rows.copy(); rows['grp'] = C.ta.group_of(rows)
    rows['loo_s'] = sur.loo_rmse; rows['loo_e'] = np.sqrt((loo_e**2).mean(0)); rows['spread'] = Y.std(0)
    print(rows.groupby('grp')[['loo_s', 'loo_e', 'spread']].median().round(3).to_string())
    out[name] = dict(p=emu.p, trunc=emu.trunc, inflation=infl, z2_median=float(np.median(z2)),
                     loo_rmse_surface=float(np.sqrt((loo_s**2).mean())), loo_rmse_emulator=float(np.sqrt((loo_e**2).mean())))
json.dump(out, open(C.OUT / 'surrogates.json', 'w'), indent=1)
