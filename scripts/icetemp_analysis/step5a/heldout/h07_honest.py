"""HONEST ERROR MODEL: misfit ratio r^T Sigma^-1 r / n at each target's best run (real model output,
Sigma = D + S, + month-to-month term and date spread for T3; no emulator term), at the five best
runs, at run 13, and at the in-sample fit through the surrogates (with E). Split into the entity
level (GLS entity means, one dimension per entity) and the within-entity rest; for T3 also the
point means vs month-to-month (target_analysis.misfit_ratio)."""
import sys, pickle, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

pd.set_option('display.width', 250, 'display.max_columns', 40)
D = H.load_grid()
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
keep = design['keep'].to_numpy()
R32 = D['runs32']['noE']
n_ent = len(D['ents'])
rec, best = [], {}
for tg in ('T1', 'T2', 'T3'):
    t = T[tg]
    tot = R32[tg]['total'].sum(1)
    lev = R32[tg]['level'].sum(1)
    pen_full = H.penalty(th3[:, :2]) + ((th3[:, 2] - 1.0) / 0.25) ** 2
    order = np.argsort(tot)
    order_post = np.argsort(tot + pen_full)
    best[tg] = int(order[0])
    for i in range(32):
        r = dict(target=tg, run=i, f_ref=th3[i, 0], f_ins=th3[i, 1], f_adv=th3[i, 2], kept=bool(keep[i]),
                 rank_chi2=int(np.flatnonzero(order == i)[0]) + 1, rank_post=int(np.flatnonzero(order_post == i)[0]) + 1,
                 n=t.n, ratio=tot[i] / t.n, entity_level=lev[i] / n_ent,
                 within_entity=(tot[i] - lev[i]) / (t.n - n_ent) if t.n > n_ent else np.nan,
                 rmse_T2=float(D['runs32']['rmse2'][i]))
        if tg == 'T3':
            mr = H.ta.misfit_ratio(t.Y[i] - t.y, t.Sigma, t.rows)
            r.update(point_means=mr['level'], month_to_month=mr['month_to_month'], check_ratio=mr['ratio'])
        rec.append(r)
R = pd.DataFrame(rec)
print('== misfit ratio at each target\'s best run (smallest misfit among the 32 real runs) ==')
for tg in ('T1', 'T2', 'T3'):
    sub = R[R.target == tg].sort_values('rank_chi2')
    print(f'\n{tg}: best run {best[tg]}; five best by misfit:')
    cols = ['run', 'f_ref', 'f_ins', 'f_adv', 'kept', 'rank_chi2', 'rank_post', 'ratio', 'entity_level', 'within_entity', 'rmse_T2']
    if tg == 'T3':
        cols += ['point_means', 'month_to_month']
    print(sub.head(5)[cols].round(3).to_string(index=False))
    r13 = sub[sub.run == 13].iloc[0]
    print(f'   run 13: ratio {r13.ratio:.3f} (rank {r13.rank_chi2}); median over the 32 runs {sub.ratio.median():.3f}; '
          f'range {sub.ratio.min():.3f}-{sub.ratio.max():.3f}; runs with ratio in 0.5-2: {((sub.ratio >= 0.5) & (sub.ratio <= 2)).sum()} of 32')
# at the in-sample fits through the surrogates (Sigma incl. E)
HO = pickle.load(open(H.OUT / 'h04_heldout.pkl', 'rb'))
INS = HO['INS']
print('\n== at the in-sample fitted parameters (surrogate; Sigma incl. the emulator term E where marked) ==')
G = D['G']
rows = []
for m, ev in (('RS', 'E'), ('RS', 'noE'), ('EMU', 'E'), ('EMU', 'noE')):
    for tg in ('T1', 'T2', 'T3'):
        i = INS[(INS.method == m) & (INS.E == ev) & (INS.target == tg)].iloc[0]
        k = int(np.argmin(np.abs(G[:, 0] - i.f_ref) + np.abs(G[:, 1] - i.f_ins)))
        p = D['methods'][m][ev][tg]
        tot, lev = p['total'][k].sum(), p['level'][k].sum()
        rows.append(dict(method=m, E=ev, target=tg, f_ref=i.f_ref, f_ins=i.f_ins, ratio=tot / T[tg].n,
                         entity_level=lev / n_ent,
                         within_entity=(tot - lev) / (T[tg].n - n_ent) if T[tg].n > n_ent else np.nan))
print(pd.DataFrame(rows).round(3).to_string(index=False))
R.to_pickle(H.OUT / 'h07_honest.pkl')
json.dump(dict(best=best, at_fit=rows), open(H.OUT / 'h07_honest.json', 'w'), indent=1, default=float)
