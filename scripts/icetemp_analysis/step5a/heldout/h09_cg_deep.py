"""Sensitivity: Colle Gnifetti below 40 m (Monte Rosa entities at band >= 4400 m, point depth > 40 m;
the accepted pre-1942 climate limitation) left out of every fit. Consistency distance and held-out
RMSE (scored on all T2 points of the left-out system, and without those points), response surface
with E."""
import sys, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

S = H.Setup(methods=('RS',))
G, pen = S.G, S.pen
t2, t3 = S.T['T2'], S.T['T3']
cg2 = ((t2.rows.system_short == 'Monte Rosa') & (t2.rows.band >= 4400) & (t2.rows.depth > 40)).to_numpy()
cg3 = ((t3.rows.system_short == 'Monte Rosa') & (t3.rows.band >= 4400) & (t3.rows.depth > 40)).to_numpy()
print('CG > 40 m: T2 points', cg2.sum(), 'T3 rows', cg3.sum(), 'entities', sorted(t2.rows[cg2].entity.unique()))
sys_of = dict(zip(S.D['ents'], S.D['sys_ent']))
systems = sorted(set(S.D['sys_ent']))
M2 = S.pred_grid('RS', 'T2')
out = {}
for lab, k2, k3 in (('all rows', np.ones(t2.n, bool), np.ones(t3.n, bool)), ('without CG > 40 m', ~cg2, ~cg3)):
    C1, e1 = S.level_t1('RS', k2)
    p2, e2 = S.parts_sub('RS', 'T2', k2)
    p3, e3 = S.parts_sub('RS', 'T3', k3)
    objs = {'T1': (C1, e1), 'T2': (p2['total'], e2), 'T3': (p3['total'], e3),
            'level': (C1, e1), 'T2 shape': (p2['total'] - p2['level'], e2), 'T3 shape': (p3['total'] - p3['level'], e3)}
    fits = {}
    for name, (C, ents) in objs.items():
        k = int(np.argmin(C.sum(1) + pen))
        th = G[k]
        if name in ('T1', 'level'):
            Fe, _ = S.level_t1('RS', k2, th=th)
            F = Fe.sum(0)
        else:
            tg = name[:2]
            fp, _ = S.parts_sub('RS', tg, k2 if tg == 'T2' else k3, th=th)
            F = H.part_of(fp, 'shape' if 'shape' in name else 'total').sum(0)
        fits[name] = dict(th=th, C=H.post_cov(F))
    d2 = H.mahal(fits['level']['th'], fits['level']['C'], fits['T2 shape']['th'], fits['T2 shape']['C'])
    d3 = H.mahal(fits['level']['th'], fits['level']['C'], fits['T3 shape']['th'], fits['T3 shape']['C'])
    # held-out
    ho = {}
    for tg in ('T1', 'T2', 'T3'):
        C, ents = objs[tg]
        sysv = np.array([sys_of[e] for e in ents])
        sse_all = sse_kept = n_all = n_kept = 0.0
        for h in systems:
            k = int(np.argmin(C @ (sysv != h).astype(float) + pen))
            te = t2.system == h
            r2 = (M2[k] - t2.y) ** 2
            sse_all += r2[te].sum(); n_all += te.sum()
            sse_kept += r2[te & ~cg2].sum(); n_kept += (te & ~cg2).sum()
        ho[tg] = dict(all_points=float(np.sqrt(sse_all / n_all)), without_cg_points=float(np.sqrt(sse_kept / n_kept)))
    out[lab] = dict(fits={k: dict(th=list(map(float, v['th'])), **H.summarize_post(v['C'])) for k, v in fits.items()},
                    d_T2=d2, d_T3=d3, heldout=ho)
    print(f'\n-- {lab} --')
    for k, v in fits.items():
        s = H.summarize_post(v['C'])
        print(f'   {k:9s} f_ref {v["th"][0]:.3f} f_ins {v["th"][1]:.3f} (post sd {s["sd_ref"]:.2f}/{s["sd_ins"]:.2f})')
    print(f'   distance level-shape: T2 {d2:.2f} sd, T3 {d3:.2f} sd')
    print('   held-out RMSE [C]:', {k: {kk: round(vv, 3) for kk, vv in v.items()} for k, v in ho.items()})
    print('   T2-T1, T3-T1 (all points):', round(ho['T2']['all_points'] - ho['T1']['all_points'], 3),
          round(ho['T3']['all_points'] - ho['T1']['all_points'], 3))
json.dump(out, open(H.OUT / 'h09_cg_deep.json', 'w'), indent=1, default=float)
