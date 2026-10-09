"""Shared measurement error for T3 rows (revise stage): the T3 misfit ratio at design 13 with D
diagonal (as used), fully shared within a point, and half shared. Checks whether Sigma stays
factorisable and how much the ratio and its month-to-month part change."""
import sys

import numpy as np

sys.path.insert(0, '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/scripts/icetemp_analysis')
import target_analysis as ta  # noqa: E402

obs = ta.load_obs()
s2, r7 = ta.load_model('alps_s2'), ta.load_model('alps_r7')
common = ta.common_rows([s2, r7], obs)
t3 = ta.target('T3r', obs, s2, common)
rows = t3.rows
V = np.diag(ta.date_spread_var(rows, obs, s2['Y3'][13]))
S = ta.s_model_cov(rows, seasonal=ta.seasonal_cov)
s = rows.sigma.to_numpy(float)
t2 = rows.t2.to_numpy(int)
same_pt = t2[:, None] == t2[None, :]
r = t3.Y[13] - t3.y
for share in (0.0, 0.5, 0.9, 1.0):
    D = (1 - share) * np.diag(s ** 2) + share * np.where(same_pt, np.outer(s, s), 0.0)
    C = D + S + V
    bad = []
    for e in np.unique(rows.e):
        j = np.flatnonzero(rows.e.to_numpy() == e)
        try:
            np.linalg.cholesky(C[np.ix_(j, j)])
        except np.linalg.LinAlgError:
            bad.append(e)
    keep = ~np.isin(rows.e.to_numpy(), bad)
    m = ta.misfit_ratio(r[keep], C[np.ix_(keep, keep)], rows[keep].reset_index(drop=True))
    ev = min(np.linalg.eigvalsh(C[np.ix_(j, j)]).min() for j in
             [np.flatnonzero(rows.e.to_numpy() == e) for e in np.unique(rows.e) if e not in bad])
    print(f'share {share:.1f}: entities not factorisable {len(bad)} {[rows.entity[rows.e == b].iloc[0] for b in bad]}; '
          f'ratio {m["ratio"]:.3f} (point means {m["level"]:.3f}, month-to-month {m["month_to_month"]:.3f}) on {m["n"]} rows; '
          f'smallest eigenvalue of a block {ev:.2e}')
