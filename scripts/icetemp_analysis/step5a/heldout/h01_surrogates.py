"""Surrogate accuracy: response-surface forms and the SVD+GP emulator, leave-one-run-out over the
29 runs with f_ref >= 0.1, measured against the error model (whitened) and on the T2 RMSE ruler."""
import sys, pickle, warnings
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

tm = H.timer()
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
keep = design['keep'].to_numpy()
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
B = {k: H.target_blocks(t) for k, t in T.items()}
print('targets and blocks', tm())
res = {}
for name in ('T2', 'T3'):
    t = T[name]
    Y = t.Y[keep]
    loos = {}
    for form in ('L3', 'Q2L', 'Q3', 'C2L', 'C2Q'):
        loos[form] = H.RespSurf(form).fit(th3[keep], Y).loo
    loos['EMU999'] = H.emu_loo(th3[keep], Y, 0.999)
    loos['EMU99'] = H.emu_loo(th3[keep], Y, 0.99)
    print(name, 'surrogates done', tm())
    for k, E in loos.items():
        p = B[name].parts(E)
        tot = p['total'].sum(1) / t.n
        lev = p['level'].sum(1)
        line = (f'{name} {k:7s} LOO rmse {np.sqrt((E**2).mean()):.3f} C (median over runs '
                f'{np.median(np.sqrt((E**2).mean(1))):.3f}, max run {np.sqrt((E**2).mean(1)).max():.3f}) | '
                f'whitened e^T S^-1 e / n: median {np.median(tot):.3f}, mean {tot.mean():.3f}, max {tot.max():.3f}')
        if 'points' in p:
            m2m = (p['total'] - p['points']).sum(1) / (t.n - len(np.unique(t.pt)))
            line += f' | month-to-month part median {np.median(m2m):.3f} max {m2m.max():.3f}'
        if name == 'T2':
            E1 = E @ T['T1'].A2.T
            p1 = B['T1'].parts(E1)
            t1w = p1['total'].sum(1) / T['T1'].n
            line += f' | T1 pooled: rmse {np.sqrt((E1**2).mean()):.3f}, whitened median {np.median(t1w):.3f} max {t1w.max():.3f}'
            # error on the RMSE ruler: RMSE(true + e) - RMSE(true)
            Ytrue = Y
            r_true = np.sqrt(((Ytrue - t.y) ** 2).mean(1))
            r_pred = np.sqrt(((Ytrue + E - t.y) ** 2).mean(1))
            line += f' | T2-RMSE error: median |d| {np.median(np.abs(r_pred - r_true)):.3f}, max {np.abs(r_pred - r_true).max():.3f}'
        print(line)
    res[name] = loos
pickle.dump(dict(loos=res, keep=keep), open(H.OUT / 'h01_loo.pkl', 'wb'))
print('done', tm())
