import sys, warnings, time
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
t1, t2, t3 = T['T1'], T['T2'], T['T3']
print(t1.rows.groupby('system_short').agg(n=('e', 'size'), ents=('entity', lambda s: ', '.join(s))).to_string())
print(t2.rows.groupby('system_short').agg(points=('e', 'size'), ents=('e', 'nunique')))
keep = design['keep'].to_numpy(); th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
e = H.SvdGp(0.999).fit(th3[keep], t3.Y[keep])
print('p T3', e.p)
G, fr, fi = H.make_grid((0.10, 0.30), (0.98, 1.79), 0.01)
t0 = time.time(); M = e.predict(G); print('emu predict', G.shape, time.time() - t0)
e2 = H.SvdGp(0.999).fit(th3[keep], t2.Y[keep]); print('p T2', e2.p)
# T1 Sigma diagonal?
off = t1.Sigma - np.diag(np.diag(t1.Sigma)); print('T1 Sigma off-diag max', np.abs(off).max())
print('design range kept', th3[keep].min(0), th3[keep].max(0))
