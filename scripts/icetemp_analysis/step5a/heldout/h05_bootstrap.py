"""90 % intervals of the held-out RMSE differences T2 - T1 and T3 - T1 (paired), three ways:
 (a) MAIN: boreholes (entities) resampled within each glacier system, folds refitted every time
     (captures both which boreholes train the fit and which are scored; boreholes are the
     independent units, their points are not);
 (b) points resampled within each held-out system, fold fits kept (the task's example);
 (c) the six folds resampled (as the 30 Sep prototype)."""
import sys, pickle, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

tm = H.timer()
D = H.load_grid()
HO = pickle.load(open(H.OUT / 'h04_heldout.pkl', 'rb'))
L = HO['L']
sys_ent, npts = D['sys_ent'], D['n_pts_ent']
systems = sorted(set(sys_ent))
n_ent = len(sys_ent)
NB = 2000
rng = np.random.default_rng(20261009)
groups = {s: np.flatnonzero(sys_ent == s) for s in systems}
W = np.zeros((NB, n_ent))
for b in range(NB):
    for s, idx in groups.items():
        np.add.at(W[b], rng.choice(idx, size=len(idx), replace=True), 1.0)

obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
t2 = T['T2']
keep = design['keep'].to_numpy()
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
sur = {'RS': H.RespSurf('C2Q').fit(th3[keep], t2.Y[keep]), 'RSQ3': H.RespSurf('Q3').fit(th3[keep], t2.Y[keep]),
       'EMU': H.SvdGp(0.999).fit(th3[keep], t2.Y[keep])}
fold_of_pt = t2.system

out = {}
for m, ev in (('RS', 'E'), ('EMU', 'E'), ('RSQ3', 'E'), ('DISC', 'noE')):
    rm = {}
    for target in ('T1', 'T2', 'T3'):
        C, pen, sse, TH = H.objective_arrays(D, m, ev, target)
        SSE = np.zeros((NB, len(systems)))
        N = np.zeros((NB, len(systems)))
        for j, h in enumerate(systems):
            te = sys_ent == h
            for c0 in range(0, NB, 250):
                Wtr = W[c0:c0 + 250] * (~te)[None, :]
                k = np.argmin(C @ Wtr.T + pen[:, None], axis=0)
                SSE[c0:c0 + 250, j] = (W[c0:c0 + 250][:, te] * sse[k][:, te]).sum(1)
            N[:, j] = W[:, te] @ npts[te]
        rm[target] = np.sqrt(SSE.sum(1) / N.sum(1))
    # (b) points within folds, fold fits fixed; (c) folds
    sub = L[(L.method == m) & (L.E == ev)]
    e2 = {}
    for target in ('T1', 'T2', 'T3'):
        v = np.zeros(t2.n)
        for h in systems:
            r = sub[(sub.target == target) & (sub.heldout == h)].iloc[0]
            pt = fold_of_pt == h
            if m == 'DISC':
                pred = t2.Y[int(r.run)]
            else:
                pred = sur[m].predict(np.array([[r.f_ref, r.f_ins]]))[0]
            v[pt] = (pred[pt] - t2.y[pt]) ** 2
        e2[target] = v
    rngb = np.random.default_rng(7)
    bp = {t: np.zeros(NB) for t in e2}
    pts_of = {h: np.flatnonzero(fold_of_pt == h) for h in systems}
    for b in range(NB):
        draw = np.concatenate([rngb.choice(i, size=len(i), replace=True) for i in pts_of.values()])
        for t in e2:
            bp[t][b] = np.sqrt(e2[t][draw].mean())
    fs = sub.pivot(index='heldout', columns='target', values='sse').loc[systems]
    fn = sub.groupby('heldout').n.first().loc[systems].to_numpy()
    bf = {t: np.zeros(NB) for t in e2}
    rngf = np.random.default_rng(11)
    for b in range(NB):
        s = rngf.choice(len(systems), size=len(systems), replace=True)
        for t in e2:
            bf[t][b] = np.sqrt(fs[t].to_numpy()[s].sum() / fn[s].sum())
    res = {}
    point = {t: float(np.sqrt(sub[sub.target == t].sse.sum() / sub[sub.target == t].n.sum())) for t in e2}
    for scheme, R in (('a_boreholes_refit', rm), ('b_points_in_folds', bp), ('c_folds', bf)):
        for d, ref in (('T2', 'T1'), ('T3', 'T1'), ('T3', 'T2')):
            diff = R[d] - R[ref]
            res[f'{scheme} {d}-{ref}'] = dict(estimate=point[d] - point[ref], p05=float(np.percentile(diff, 5)),
                                          p50=float(np.percentile(diff, 50)), p95=float(np.percentile(diff, 95)),
                                          frac_above_0p10=float((diff > 0.10).mean()), frac_below_0=float((diff < 0).mean()))
        for t in ('T1', 'T2', 'T3'):
            res[f'{scheme} {t}'] = dict(estimate=point[t], p05=float(np.percentile(R[t], 5)), p95=float(np.percentile(R[t], 95)))
    out[f'{m}/{ev}'] = res
    print(f'\n== {H.METHOD_LABEL[(m, ev)]} ({tm()}) ==')
    print(pd.DataFrame(res).T.round(3).to_string())
json.dump(out, open(H.OUT / 'h05_bootstrap.json', 'w'), indent=1)
