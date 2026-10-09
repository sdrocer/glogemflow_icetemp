"""REML fits of candidate month-to-month kernels at design 13 (alps_s2), T3 month rows."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import sys, time, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import multiprocessing as mp
import numpy as np, pandas as pd
import target_analysis as ta
import s2lib as L
import s2fit as F
OUT = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
DESIGN = int(sys.argv[1]) if len(sys.argv) > 1 else 13
r = Y[DESIGN] - y
B = F.Blocks(rows)
CANDS = {
  'M1 exp only':              dict(a=1.0, ze=5.0, l0=2.0, tr=2.0),
  'M2 exp, l grows':          dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, tr=2.0),
  'M3 periodic+exp':          dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tr=2.0),
  'M4 quasi-periodic+exp':    dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tp=36.0, tr=2.0),
  'M5 M3 + slow drift':       dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tr=2.0, aC=0.5, tC=120.0, lC=50.0),
  'M6 M4 + slow drift':       dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tp=36.0, tr=2.0, aC=0.5, tC=120.0, lC=50.0),
  'M7 M5, drift fades in z':  dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tr=2.0, aC=0.5, zC=20.0, tC=120.0, lC=50.0),
  'M8 M6, drift fades in z':  dict(a=1.0, ze=5.0, l0=1.0, l1=0.3, w=0.5, tp=36.0, tr=2.0, aC=0.5, zC=20.0, tC=120.0, lC=50.0),
}


def job(name):
    t0 = time.time()
    p, ll, k = F.fit(B, r, CANDS[name])
    return name, p, ll, k, time.time() - t0


if __name__ == '__main__':
    print(f'design {DESIGN}: blocks {len(B.blocks)}, rows {B.n}, points {B.p}, contrasts {B.n - B.p}', flush=True)
    with mp.get_context('fork').Pool(len(CANDS)) as pool:
        res = pool.map(job, list(CANDS))
    out = {}
    for name, p, ll, k, s in res:
        out[name] = (p, ll, k)
        print(f'{name:26s} REML {ll:9.2f} k {k:2d} AIC {-2*ll+2*k:9.2f} ({s:.0f} s) ' +
              ' '.join(f'{a}={b:.4g}' for a, b in p.items()), flush=True)
    pickle.dump(out, open(OUT + f'reml_d{DESIGN}.pkl', 'wb'))
