import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import s2g_within as G
P, R = G.P, G.R
groups = {'all': np.ones(len(P), bool), 'Corvatsch': (P.system_short == 'Corvatsch').to_numpy(),
          'Monte Rosa': (P.system_short == 'Monte Rosa').to_numpy(),
          'all but Monte Rosa': (P.system_short != 'Monte Rosa').to_numpy(),
          'all but GG94 (Grenz 4105-4305)': ~P.entity.isin(['Grenzgletscher@4105', 'Grenzgletscher@4205', 'Grenzgletscher@4305']).to_numpy()}
full = P.copy()
for name, m in groups.items():
    G.P = full[m].reset_index(drop=True)
    p, cl, _ = G.fit(R[13], 'exp')
    sub = G.P
    print(f'{name:32s} pairs {len(sub):5d} pts {sub.t2.nunique():4d} ' + ' '.join(f'{k}={v:.3g}' for k, v in p.items())
          + f'  v(1m)={p["a"]*np.exp(-1/p["zs"]):.2f} v(5m)={p["a"]*np.exp(-5/p["zs"]):.2f} v(10m)={p["a"]*np.exp(-10/p["zs"]):.3f}')
G.P = full
