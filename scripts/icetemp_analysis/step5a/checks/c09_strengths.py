"""Supporting numbers: how strongly f_ref and f_ins act on the row sets the checks move (surface
Jacobian at the baseline optimum, C per prior sd: f_ref 0.2598, f_ins 0.3259), and, for (c), the
chi2 the one-sided treatment forgives at the optimum. Writes c09_strengths.log."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
mel = pickle.load(open(C.OUT / 'c05_melting.pkl', 'rb'))['near']
sn = pickle.load(open(C.OUT / 'c04c_snow_fit.pkl', 'rb'))
HEF, CG = 'Hintereisferner tongue, ice (< 3000 m)', 'Colle Gnifetti (>= 4400 m)'
for name in ('T2', 'T3r'):
    c = ctx[name]
    s = base[name]['surface']
    th0 = np.array([s['map_f_ref'], s['map_f_ins'], 1.0])
    lin = C.linearised(c.sur, c.y, c.Sigma, th0)
    J = lin['J']
    R = c.rows
    sets = {'all': np.ones(len(R), bool), 'top 3 m': (R.depth <= 3).to_numpy(), '3-15 m': ((R.depth > 3) & (R.depth <= 15)).to_numpy(),
            'below 15 m': (R.depth > 15).to_numpy(), 'HEF tongue': (R.grp == HEF).to_numpy(),
            'CG > 40 m': ((R.grp == CG) & (R.depth > 40)).to_numpy(), 'CG <= 40 m': ((R.grp == CG) & (R.depth <= 40)).to_numpy(),
            'at melting point': mel[name], 'old (before 1979)': (R.ym_last < 197901).to_numpy(),
            'shallow under snow': sn['p2_snow'] if name == 'T2' else sn['r3_snow_shallow'],
            'label firn': (R.firn_label == 'firn').to_numpy(), 'label ice': (R.firn_label == 'ice').to_numpy()}
    rec = []
    for k, m in sets.items():
        Jm = J[m]
        cos = float(Jm[:, 0] @ Jm[:, 1] / (np.linalg.norm(Jm[:, 0]) * np.linalg.norm(Jm[:, 1]) + 1e-300))
        rec.append(dict(rows=k, n=int(m.sum()), rms_J_f_ref=float(np.sqrt((Jm[:, 0] ** 2).mean())),
                        rms_J_f_ins=float(np.sqrt((Jm[:, 1] ** 2).mean())), cos_ref_ins=cos))
    print(f'== {name}: sensitivity at the surface optimum ({th0[0]:.3f}, {th0[1]:.3f}), C per prior sd ==')
    print(pd.DataFrame(rec).round(3).to_string(index=False))
    # (c) chi2 forgiven by the one-sided treatment at the optimum
    G = c.sur.predict(th0)[0]
    r = G - c.y
    r1 = C.hinge(r[None, :], c.y, mel[name], c.cap)[0]
    q0 = C.chi2_fixed(r[None, :], c.Sigma)[0][0]; q1 = C.chi2_fixed(r1[None, :], c.Sigma)[0][0]
    zero = mel[name] & (r1 == 0) & (r > 0)
    capd = mel[name] & (r1 != r) & ~zero
    print(f'(c) {name}: chi2 at the optimum ordinary {q0:.2f}, one-sided {q1:.2f} (change {q1 - q0:+.2f}); of the {mel[name].sum()} rows at the '
          f'melting point, {zero.sum()} have the model at or above the reading (residual +{r[zero].min():.3f} to +{r[zero].max():.3f} C, '
          f'set to 0) and {capd.sum()} have a reading above PMP and a colder model (residual measured from PMP instead: '
          f'{(r1 - r)[capd].min():.3f} to {(r1 - r)[capd].max():.3f} C smaller); the rest unchanged\n')
