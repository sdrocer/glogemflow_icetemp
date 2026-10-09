"""Held-out RMSE by glacier system (6 folds): fit (f_ref, f_ins) on five systems by the target's GLS
misfit + prior penalty, score the left-out system's T2 points at the fitted parameters; pooled over
folds. Also the in-sample (all six systems) fit of each target. Saves h04_heldout.pkl."""
import sys, pickle, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

pd.set_option('display.width', 250, 'display.max_columns', 40)
D = H.load_grid()
sys_ent, npts = D['sys_ent'], D['n_pts_ent']
systems = sorted(set(sys_ent))
N_all = npts.sum()
recs, ins = [], []
for m, ev in H.METHODS + (('RS', 'E-rmse'), ('EMU', 'E-rmse'), ('DISC', 'noE-rmse')):
    rmse_fit = ev.endswith('-rmse')
    ev0 = ev.replace('-rmse', '')
    for target in (('T2',) if rmse_fit else ('T1', 'T2', 'T3')):
        C, pen, sse, TH = H.objective_arrays(D, m, ev0, target)
        if rmse_fit:                     # reference: plain sum of squares on T2 points, no prior
            C, pen = sse, np.zeros(len(TH))
        tname = 'RMSE fit' if rmse_fit else target
        k, J = H.argmin_w(C, pen, np.ones(C.shape[1]))
        ins.append(dict(method=m, E=ev0, target=tname, f_ref=TH[k, 0], f_ins=TH[k, 1], edge=H.on_edge(TH[k]),
                        rmse_T2=np.sqrt(sse[k].sum() / N_all), chi2=float(C[k].sum()), run=int(D['discrete']['runs'][k]) if m == 'DISC' else -1))
        for h in systems:
            w = (sys_ent != h).astype(float)
            k, J = H.argmin_w(C, pen, w)
            te = sys_ent == h
            recs.append(dict(method=m, E=ev0, target=tname, heldout=h, f_ref=TH[k, 0], f_ins=TH[k, 1],
                             edge=H.on_edge(TH[k]), k=k, sse=float(sse[k, te].sum()), n=int(npts[te].sum()),
                             run=int(D['discrete']['runs'][k]) if m == 'DISC' else -1))
L = pd.DataFrame(recs)
L['rmse'] = np.sqrt(L.sse / L.n)
INS = pd.DataFrame(ins)
pool = L.groupby(['method', 'E', 'target']).agg(sse=('sse', 'sum'), n=('n', 'sum')).reset_index()
pool['heldout_rmse'] = np.sqrt(pool.sse / pool.n)
print('== pooled held-out RMSE on the T2 points of the left-out system [C] ==')
piv = pool.pivot_table(index=['method', 'E'], columns='target', values='heldout_rmse')
for t in ('T2', 'T3'):
    piv[f'{t}-T1'] = piv[t] - piv['T1']
print(piv.round(3).to_string())
print('\n== in-sample fit (all six systems): parameters and T2 RMSE (771 points) ==')
print(INS.round(3).to_string())
for m, ev in (('RS', 'E'), ('EMU', 'E'), ('DISC', 'noE')):
    sub = L[(L.method == m) & (L.E == ev)]
    print(f'\n== folds, {H.METHOD_LABEL[(m, ev)]}: held-out RMSE [C] and fitted f_ref / f_ins ==')
    a = sub.pivot(index='heldout', columns='target', values='rmse').round(2)
    b = sub.assign(th=sub.f_ref.round(2).astype(str) + ' / ' + sub.f_ins.round(2).astype(str)
                   + np.where(sub.edge, ' (edge)', '')).pivot(index='heldout', columns='target', values='th')
    a['points'] = sub.groupby('heldout').n.first()
    print(a.to_string())
    print(b.to_string())
    if m == 'DISC':
        print(sub.pivot(index='heldout', columns='target', values='run').to_string())
pickle.dump(dict(L=L, INS=INS, pool=pool), open(H.OUT / 'h04_heldout.pkl', 'wb'))
