"""Probe: Emulator (SVD 0.999 + PCGPwM) on target rows, timing, predictive variance, LOO inflation."""
import sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
import target_analysis as ta
from icetemp.calibration import emulator as EM

EM.Emulator.assemble_matrix = lambda self, outs: np.column_stack([np.asarray(o, float) for o in outs])

obs = ta.load_obs(); m = ta.load_model('alps_s2'); r7 = ta.load_model('alps_r7')
c = ta.common_rows([m, r7], obs)
d = ta.load_design('alps_s2'); th = d[['f_ref', 'f_ins', 'f_adv']].to_numpy()
for name in ('T2', 'T3r'):
    t = ta.target(name, obs, m, c)
    for keep in (np.ones(32, bool), d.keep.to_numpy()):
        t0 = time.time()
        emu = EM.Emulator(glaciers=[], explained_variance=0.999)
        emu.fit(th[keep], list(t.Y[keep]))
        t1 = time.time() - t0
        mu, L = emu.predict_factor(th[13])
        mu2, var2 = emu.predict(th[13])
        th_off = th[13].copy(); th_off[2] = 1.0
        _, Loff = emu.predict_factor(th_off)
        mid = 0.5 * (th[13] + th[14])
        _, Lmid = emu.predict_factor(mid)
        print(f'{name} n={keep.sum()} p={emu.p} fit {t1:.1f}s trunc {emu.truncation_variance():.2e} '
              f'| at d13: max|mu-Y| {np.abs(mu - t.Y[13]).max():.2e} trace(LL^T)/n {np.sum(L**2)/len(mu):.2e} '
              f'| f_adv=1: {np.sum(Loff**2)/len(mu):.2e} | mid13-14: {np.sum(Lmid**2)/len(mu):.2e}')
        if keep.all():
            t0 = time.time()
            emu.calibrate_variance(th[keep], list(t.Y[keep]))
            print(f'   calibrate_variance {time.time()-t0:.0f}s: loo_z_std {emu.loo_z_std_:.2f} infl flat {emu.variance_inflation_:.1f} '
                  f'curve thresholds X {np.round(emu._variance_inflation_curve.X_thresholds_, 3)} y {np.round(emu._variance_inflation_curve.y_thresholds_, 1)} '
                  f'trust {emu.max_trusted_distance_:.3f} spacing {emu.design_spacing_:.3f}')
            _, Lc = emu.predict_factor(th_off)
            print(f'   after calibration: f_adv=1 trace/n {np.sum(Lc**2)/len(mu):.2e}  nn dist {emu._nn_distance(th_off)}')
