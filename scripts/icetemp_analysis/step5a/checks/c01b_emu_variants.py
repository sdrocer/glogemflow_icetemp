"""Why the emulator posterior runs to f_ins = 2: E term and extrapolation. Variants at baseline."""
import sys, time
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ctx = C.context()
for name, c in ctx.items():
    TH = C.grid()
    m = np.ones(len(c.y), bool)
    emu = c.emu
    G, cvar = emu.predict(TH)
    R = G - c.y[None, :]
    c2_mean, ld0 = C.chi2_fixed(R, c.Sigma)
    c2_E, ldE = C.chi2_emu(R, c.Sigma, emu.basis, cvar, emu.trunc)
    # design box and trust radius (pipeline rule) in 3D standardized
    th_tr = c.d.theta[c.d.keep]
    sd = th_tr.std(axis=0)
    Z = th_tr / sd
    pair = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1)); np.fill_diagonal(pair, np.inf)
    spacing = np.median(pair.min(axis=1)); maxd = 2 * spacing
    nn = np.sqrt((((TH / sd)[:, None, :] - Z[None, :, :]) ** 2).sum(-1)).min(axis=1)
    lo, hi = th_tr.min(0), th_tr.max(0)
    inbox = (TH[:, 0] >= lo[0]) & (TH[:, 0] <= hi[0]) & (TH[:, 1] >= lo[1]) & (TH[:, 1] <= hi[1])
    trusted = inbox & (nn <= maxd)
    print(f'{name}: design box f_ref {lo[0]:.3f}-{hi[0]:.3f} f_ins {lo[1]:.3f}-{hi[1]:.3f}; spacing {spacing:.3f}, trust radius {maxd:.3f}; '
          f'grid trusted {trusted.mean():.2%}')
    lp_prior = C.log_prior(TH)
    for lab, lp, mask in (('mean only, full grid', -0.5 * c2_mean + lp_prior, np.ones(len(TH), bool)),
                          ('mean only, trusted', -0.5 * c2_mean + lp_prior, trusted),
                          ('with E (infl %.2f), full grid' % emu.infl, -0.5 * (c2_E + ldE) + lp_prior, np.ones(len(TH), bool)),
                          ('with E, trusted', -0.5 * (c2_E + ldE) + lp_prior, trusted)):
        lpm = np.where(mask, lp, -np.inf)
        s = C.summarize_grid(TH, lpm)
        k = int(np.argmax(lpm))
        print(f'   {lab:32s}: MAP {TH[k,0]:.2f}/{TH[k,1]:.2f} mean {s["mean_f_ref"]:.3f}/{s["mean_f_ins"]:.3f} corr {s["corr"]:+.2f} '
              f'spread {s["spread_f_ref"]:.2f}/{s["spread_f_ins"]:.2f}; chi2 mean-only at MAP {c2_mean[k]:.1f}, with E {c2_E[k]:.1f}, '
              f'logdet-ld0 {ldE[k]-ld0:.1f}')
    # E size: mean coefficient-variance-induced row variance at the surface MAP vs at f_ins 2
    for th in ([0.33, 0.76, 1.0], [0.43, 2.0, 1.0], [0.52, 1.09, 1.0]):
        _, cv = emu.predict(np.array(th))
        v = (emu.basis ** 2) @ cv[0]
        print(f'   E at {th}: median row var {np.median(v):.3f}, max {v.max():.2f}, coefficient var {np.round(cv[0], 2)}')
