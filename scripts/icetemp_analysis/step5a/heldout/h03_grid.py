"""Per-entity GLS misfits (total, GLS level, GLS point means) and held-out SSE on the T2 ruler,
on a parameter grid (f_adv = 1) for each surrogate, with and without the emulator term E; and at
the 29 real runs (discrete). Saved to h03_grid.pkl."""
import sys, pickle, warnings
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

tm = H.timer()
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
t1, t2, t3 = T['T1'], T['T2'], T['T3']
keep = design['keep'].to_numpy()
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
G, fr, fi = H.make_grid((0.10, 0.30), (0.98, 1.79), 0.005)
print('grid', G.shape, fr[[0, -1]], fi[[0, -1]])
loo = pickle.load(open(H.OUT / 'h01_loo.pkl', 'rb'))['loos']
ents = np.unique(t2.ent)
assert np.array_equal(ents, np.unique(t1.ent)) and np.array_equal(ents, np.unique(t3.ent))
assert np.array_equal(t1.ent, ents)
sys_ent = t1.rows['system_short'].to_numpy(object)
n_pts_ent = np.array([(t2.ent == e).sum() for e in ents])


def sse_ent(M2):
    return np.column_stack([((M2[:, t2.ent == e] - t2.y[t2.ent == e]) ** 2).sum(1) for e in ents])


def all_parts(M1, M2, M3, E1, E2, E3):
    B1 = H.Blocks(t1.Sigma + E1, t1.ent)
    B2 = H.Blocks(t2.Sigma + E2, t2.ent)
    B3 = H.Blocks(t3.Sigma + E3, t3.ent, pts=t3.pt)
    p1 = B1.parts(M1 - t1.y)
    p2 = B2.parts(M2 - t2.y)
    p3 = B3.parts(M3 - t3.y)
    return dict(T1=p1, T2=p2, T3=p3)


out = dict(G=G, fr=fr, fi=fi, ents=ents, sys_ent=sys_ent, n_pts_ent=n_pts_ent, methods={})
zero = dict(T1=np.zeros((t1.n, t1.n)), T2=np.zeros((t2.n, t2.n)), T3=np.zeros((t3.n, t3.n)))
for m in ('RS', 'EMU', 'RSQ3'):
    if m in ('RS', 'RSQ3'):
        form = 'C2Q' if m == 'RS' else 'Q3'
        s2_ = H.RespSurf(form).fit(th3[keep], t2.Y[keep])
        s3_ = H.RespSurf(form).fit(th3[keep], t3.Y[keep])
        L2, L3 = s2_.loo, s3_.loo
    else:
        s2_ = H.SvdGp(0.999).fit(th3[keep], t2.Y[keep])
        s3_ = H.SvdGp(0.999).fit(th3[keep], t3.Y[keep])
        L2, L3 = loo['T2']['EMU999'], loo['T3']['EMU999']
    M2 = s2_.predict(G)
    M3 = s3_.predict(G)
    M1 = M2 @ t1.A2.T
    E2f = L2.T @ L2 / len(L2)
    E3f = L3.T @ L3 / len(L3)
    E = dict(T1=np.diag(np.diag(t1.A2 @ E2f @ t1.A2.T)), T2=H.block_part(E2f, t2.ent),
             T3=H.block_part(E3f, t3.ent))
    res = dict(E=E)
    for ev, EE in (('E', E), ('noE', zero)):
        res[ev] = all_parts(M1, M2, M3, EE['T1'], EE['T2'], EE['T3'])
        print(m, ev, 'parts done', tm(), flush=True)
    res['sse'] = sse_ent(M2)
    # T2 RMSE over all points and the pooled point means of the T3 surrogate (consistency check)
    res['rmse2'] = np.sqrt(((M2 - t2.y) ** 2).mean(1))
    W3 = np.zeros((t2.n, t3.n))
    cnt = t3.rows['count'].to_numpy(float)
    for p in range(t2.n):
        j = np.flatnonzero(t3.pt == p)
        W3[p, j] = cnt[j] / cnt[j].sum()
    d23 = np.abs(M3[::97] @ W3.T - M2[::97]).max()
    print(m, 'max |pooled T3 surrogate - T2 surrogate| on a grid subsample:', d23)
    res['d23'] = float(d23)
    out['methods'][m] = res
    del M1, M2, M3
# discrete: the 29 real runs
idx = np.flatnonzero(keep)
Y1 = t1.Y[idx]; Y2 = t2.Y[idx]; Y3 = t3.Y[idx]
out['discrete'] = dict(runs=idx, theta=th3[idx], noE=all_parts(Y1, Y2, Y3, zero['T1'], zero['T2'], zero['T3']),
                       sse=sse_ent(Y2), rmse2=np.sqrt(((Y2 - t2.y) ** 2).mean(1)))
# all 32 runs (for the best-run misfit ratios)
out['runs32'] = dict(theta=th3, noE=all_parts(t1.Y, t2.Y, t3.Y, zero['T1'], zero['T2'], zero['T3']),
                     rmse2=np.sqrt(((t2.Y - t2.y) ** 2).mean(1)))
pickle.dump(out, open(H.OUT / 'h03_grid.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
print('saved', tm())
