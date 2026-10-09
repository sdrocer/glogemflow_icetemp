"""Partial-point row drops in T3r: the error covariance as a submatrix of the full one (used) vs
rebuilt on the kept rows (seasonal term re-centred on the kept rows). Cases (c), (d), (g)."""
import sys, pickle, time
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np
import c3lib as C
ctx = C.context(('T3r',))
c = ctx['T3r']
mel = pickle.load(open(C.OUT / 'c05_melting.pkl', 'rb'))['near']['T3r']
sn = pickle.load(open(C.OUT / 'c04c_snow_fit.pkl', 'rb'))['r3_snow_shallow']
old = (c.rows.ym_last < 197901).to_numpy()
for lab, keep in (('(c) melting-point rows left out', ~mel), ('(d) shallow rows under snow left out', ~sn), ('(g) old rows left out', ~old)):
    t0 = time.time()
    rows = c.rows[keep].reset_index(drop=True)
    S_new = C.sigma_T3r(rows, c.d.obs, c.d.mod['Y3'][13])
    a = C.fit_rows(c, keep)
    b = C.fit_rows(c, keep, Sigma=S_new)
    for k in ('surface', 'emulator'):
        print(f'{lab:38s} {k:8s}: submatrix f_ref {a[k]["mean_f_ref"]:.3f} f_ins {a[k]["mean_f_ins"]:.3f} corr {a[k]["corr"]:+.2f} spread {a[k]["spread_f_ins"]:.2f} | '
              f'rebuilt f_ref {b[k]["mean_f_ref"]:.3f} f_ins {b[k]["mean_f_ins"]:.3f} corr {b[k]["corr"]:+.2f} spread {b[k]["spread_f_ins"]:.2f} [{time.time()-t0:.0f} s]')
