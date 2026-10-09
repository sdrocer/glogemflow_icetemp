"""Rounded proposal (from the median over the 29 runs with f_ref >= 0.1) against the design-13 fit."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2kern2 as K2
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
d = pickle.load(open(W + 'final_kernel.pkl', 'rb'))
LK = [1.2, 4.3, 200.0, 200.0]
ROUNDED = dict(a=3.0, zs=3.5, w=0.27, tau=7.0, q=0.25, hf=5.4, u0=0.45, beta=0.035, ld=50.0, l0=0, l1=0,
               lb0=0, lb1=0, lk=LK, lkb=LK, nu=1.5, nub=1.5, centre=True)
d['rounded'] = ROUNDED
pickle.dump(d, open(W + 'final_kernel.pkl', 'wb'))
for lab in ('final_mono', 'rounded'):
    p = d[lab]
    for name in ('T2', 'T3r', 'T3m', 'T3'):
        t = ta.target('T3' if name == 'T3m' else name, obs, mod)
        rr, Yt, yt = t.rows, t.Y, t.y
        if name == 'T3m':
            m = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[m].reset_index(drop=True), Yt[:, m], yt[m]
        K = K2.seasonal_matrix(rr, p, obs=obs)
        Sig = ta.error_cov(rr) + K
        extra = ''
        if name == 'T3r':
            mm0 = K2.misfit(Yt[13] - yt, Sig, rr)
            extra = f'  (without date spread {mm0["ratio"]:.3f})'
            Sig = Sig + np.diag(K2.date_spread(rr, obs, mod['Y3'][13]))
        mm = K2.misfit(Yt[13] - yt, Sig, rr)
        print(f"{lab:11s} {name:4s} ratio {mm['ratio']:.3f} level {mm['level']:.3f} contrast {mm['contrast']:.3f} "
              f"min eig K {np.linalg.eigvalsh(K).min():.2e}{extra}")
