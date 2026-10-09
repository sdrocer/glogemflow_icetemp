"""Misfit ratios with D + S only (Stage 1 smoke test reproduced) and the level/contrast split."""
import sys; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2kern as K
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
for name in ('T2', 'T3', 'T3r'):
    t = ta.target(name, obs, mod)
    S = ta.target_error_cov(t, obs)
    m = K.misfit(t.Y[13] - t.y, S, t.rows)
    print(f"{name}: n {m['n']} points {m['points']} ratio {m['ratio']:.3f} level {m['level']:.3f} contrast {m['contrast']:.3f}")
t1 = ta.target('T1', obs, mod)
S1 = ta.target_error_cov(t1, obs); r1 = t1.Y[13] - t1.y
print('T1 ratio', float(r1 @ np.linalg.solve(S1, r1)) / len(r1))
