"""Check: the shape misfit used in the fits (total minus GLS entity-level part) equals the misfit of
'rows minus their depth-weighted entity mean' with the matching covariance P Sigma P^T
(pseudo-inverse), entity by entity, for T2 and T3 at design 13 and at one grid point."""
import sys, warnings
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np
import hlib as H

obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
worst = 0.0
for tg in ('T2', 'T3'):
    t = T[tg]
    B = H.Blocks(t.Sigma, t.ent)
    r = t.Y[13] - t.y
    p = B.parts(r[None, :])
    shape_gls = (p['total'] - p['level'])[0]
    w = t.rows['weight'].to_numpy(float) * (t.rows['count'].to_numpy(float) if tg == 'T3' else 1.0)
    vals = []
    for k, i in enumerate(B.idx):
        a = w[i] / w[i].sum()
        P = np.eye(len(i)) - np.outer(np.ones(len(i)), a)
        d = P @ r[i]
        C = P @ t.Sigma[np.ix_(i, i)] @ P.T
        q = d @ np.linalg.pinv(C, rcond=1e-10, hermitian=True) @ d if len(i) > 1 else 0.0
        vals.append(q)
    diff = np.abs(np.array(vals) - shape_gls)
    rel = diff / np.maximum(1.0, shape_gls)
    worst = max(worst, rel.max())
    print(f'{tg}: shape misfit at design 13, GLS split {shape_gls.sum():.4f}, deviations from the weighted mean with '
          f'P Sigma P^T {np.sum(vals):.4f}; largest per-entity difference {diff.max():.2e}')
print('largest relative difference', worst)
