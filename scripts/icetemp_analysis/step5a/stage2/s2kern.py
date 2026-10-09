"""General month-to-month kernel for T2 / T3 / T3r rows and misfit ratios with D + S + K."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve
import target_analysis as ta

W12 = 2 * np.pi / 12.0


def ym_to_t(ym, frac):
    ym = np.asarray(ym, dtype=int)
    return (ym // 100) * 12 + ym % 100 - 1 + np.asarray(frac, dtype=float)


def row_times(rows, obs=None):
    """Per row: (times [months since year 0], weights summing to 1). T3: its month; T3r: its months
    (equal weights); T2: the point's months weighted by their reading counts (from obs['T3'])."""
    if 'kind' in rows:                                     # T3r
        return [(ym_to_t(m, f), np.full(len(m), 1.0 / len(m))) for m, f in zip(rows.months, rows.fracs)]
    if 'ym' in rows:                                       # T3
        t = ym_to_t(rows.ym.to_numpy(), rows.frac.to_numpy())
        return [(np.array([x]), np.ones(1)) for x in t]
    T3 = obs['T3']                                         # T2
    g = T3.groupby('t2')
    out = []
    for t2 in rows.t2.to_numpy():
        s = T3.iloc[g.indices[t2]]
        w = s['count'].to_numpy(float)
        out.append((ym_to_t(s.ym.to_numpy(), s.frac.to_numpy()), w / w.sum()))
    return out


def paciorek(z, l0, l1=0.0):
    l = l0 + l1 * z
    L2 = 0.5 * (l[:, None] ** 2 + l[None, :] ** 2)
    return np.sqrt(np.outer(l, l) / L2) * np.exp(-np.abs(z[:, None] - z[None, :]) / np.sqrt(L2))


def time_corr(dt, comp):
    kind = comp['kind']
    if kind == 'cos':
        return np.cos(W12 * dt) * (np.exp(-dt / comp['tp']) if comp.get('tp') else 1.0)
    if kind == 'exp':
        return np.exp(-dt / comp['tau'])
    if kind == 'gauss':
        return np.exp(-(dt / comp['tau']) ** 2)
    if kind == 'const':
        return np.ones_like(dt)
    raise ValueError(kind)


def var_z(z, comp):
    v = comp['a'] * np.exp(-z / comp['ze']) if np.isfinite(comp.get('ze', np.inf)) else np.full(len(z), comp['a'])
    return v + comp.get('b', 0.0)


def kernel_block(z, times, comps):
    """Kernel between rows of one entity: sum over components of
    sqrt(v(z_i) v(z_j)) * depth corr(z_i, z_j) * mean over the rows' times of time corr."""
    n = len(z)
    nat = [len(t) for t, _ in times]
    T = np.concatenate([t for t, _ in times])
    Wm = np.zeros((n, len(T)))
    k = 0
    for i, (t, w) in enumerate(times):
        Wm[i, k:k + len(t)] = w
        k += len(t)
    dt = np.abs(T[:, None] - T[None, :])
    K = np.zeros((n, n))
    for c in comps:
        v = var_z(z, c)
        rho = Wm @ time_corr(dt, c) @ Wm.T
        K += np.sqrt(np.outer(v, v)) * paciorek(z, c['l0'], c.get('l1', 0.0)) * rho
    return K


def seasonal_matrix(rows, comps, obs=None, times=None):
    """Full (rows x rows) kernel, zero across entities."""
    times = row_times(rows, obs) if times is None else times
    e = rows['e'].to_numpy(int); z = rows['depth'].to_numpy(float)
    K = np.zeros((len(rows), len(rows)))
    for ent in np.unique(e):
        idx = np.flatnonzero(e == ent)
        K[np.ix_(idx, idx)] = kernel_block(z[idx], [times[i] for i in idx], comps)
    return K


def misfit(r, Sigma, rows, by_point=True):
    """r^T Sigma^-1 r / n with Sigma block-diagonal by entity; also the level / contrast split
    (GLS point means: level part has p dims, contrast part n - p) and per-row contributions
    (r_i * (Sigma^-1 r)_i)."""
    e = rows['e'].to_numpy(int); t2 = rows['t2'].to_numpy(int)
    tot = lvl = 0.0; n_pts = 0
    contrib = np.zeros(len(r))
    for ent in np.unique(e):
        idx = np.flatnonzero(e == ent)
        c = cho_factor(Sigma[np.ix_(idx, idx)], lower=True, check_finite=False)
        rr = r[idx]
        Sir = cho_solve(c, rr, check_finite=False)
        contrib[idx] = rr * Sir
        tot += rr @ Sir
        pts, inv = np.unique(t2[idx], return_inverse=True)
        n_pts += len(pts)
        Z = np.zeros((len(idx), len(pts))); Z[np.arange(len(idx)), inv] = 1.0
        A = Z.T @ cho_solve(c, Z, check_finite=False)
        u = Z.T @ Sir
        lvl += u @ np.linalg.solve(A, u)
    n = len(r)
    return dict(n=n, points=n_pts, ratio=tot / n, chi2=tot, level=lvl / n_pts,
                contrast=(tot - lvl) / max(n - n_pts, 1), contrib=contrib)
