"""Full month-to-month kernel for T2 / T3 / T3r rows:
  within-borehole part  sqrt(m_i m_j v(z_i) v(z_j)) c_s(z_i, z_j) rho(dt), rho = w cos(2 pi dt/12) + (1-w) exp(-dt/tau)
                        v(z) = a exp(-z/zs); m = 1 firn/unknown, q ice; rows spanning several months
                        average rho over their months (weights)
  between-borehole part sqrt(u(z_i) u(z_j)) c_b(z_i, z_j) (b_i - bbar_p(i)) . (b_j - bbar_p(j)),
                        b = the row's borehole fractions, bbar_p their count-weighted mean over the
                        point's rows (so a point's count-weighted mean is untouched); u(z) = u0 exp(-z/zu)
  depth correlations    non-stationary Matern (nu 1/2 or 3/2) with length l(z) = l0 + l1 z (Paciorek)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve

W12 = 2 * np.pi / 12.0


def ym_to_t(ym, frac):
    ym = np.asarray(ym, dtype=int)
    return (ym // 100) * 12 + ym % 100 - 1 + np.asarray(frac, dtype=float)


def row_times(rows, obs=None):
    if 'kind' in rows:
        return [(ym_to_t(m, f), np.full(len(m), 1.0 / len(m))) for m, f in zip(rows.months, rows.fracs)]
    if 'ym' in rows:
        t = ym_to_t(rows.ym.to_numpy(), rows.frac.to_numpy())
        return [(np.array([x]), np.ones(1)) for x in t]
    T3 = obs['T3']
    g = T3.groupby('t2')
    out = []
    for t2 in rows.t2.to_numpy():
        s = T3.iloc[g.indices[t2]]
        w = s['count'].to_numpy(float)
        out.append((ym_to_t(s.ym.to_numpy(), s.frac.to_numpy()), w / w.sum()))
    return out


def matern(d, nu):
    if nu == 0.5:
        return np.exp(-d)
    if nu == 1.5:
        s = np.sqrt(3.0) * d
        return (1 + s) * np.exp(-s)
    if nu == 2.5:
        s = np.sqrt(5.0) * d
        return (1 + s + s * s / 3) * np.exp(-s)
    raise ValueError(nu)


LKNOTS = np.array([2.5, 9.0, 20.0, 40.0])


def length_z(z, l0, l1=0.0, zl=None, lk=None):
    """Depth-correlation length [m] at depth z: knots (log-linear, constant beyond), exponential
    growth l0 exp(z/zl), or linear l0 + l1 z."""
    if lk is not None:
        return np.exp(np.interp(z, LKNOTS, np.log(lk)))
    return l0 * np.exp(z / zl) if zl else l0 + l1 * z


def depth_corr(z, l0, l1, nu, zl=None, lk=None):
    l = length_z(z, l0, l1, zl, lk)
    L2 = 0.5 * (l[:, None] ** 2 + l[None, :] ** 2)
    return np.sqrt(np.outer(l, l) / L2) * matern(np.abs(z[:, None] - z[None, :]) / np.sqrt(L2), nu)


def borehole_vectors(rows):
    labs = [str(b).split('/') for b in rows.bh_label]
    names = sorted({x for l in labs for x in l})
    idx = {n: k for k, n in enumerate(names)}
    B = np.zeros((len(rows), len(names)))
    for i, l in enumerate(labs):
        for x in l:
            B[i, idx[x]] += 1.0 / len(l)
    return B


def kernel_block(sub, times, p):
    """Kernel between the rows of one entity (DataFrame sub, times list aligned)."""
    z = sub['depth'].to_numpy(float)
    n = len(z)
    ice = (sub['firn_label'].to_numpy(object) == 'ice')
    m = np.where(ice, p['q'], 1.0)
    if p.get('hf', 1.0) != 1.0:
        m = np.where(~ice & (sub['band'].to_numpy(float) >= p.get('hf_band', 4000.0)), p['hf'], m)
    v = p['a'] * np.exp(-z / p['zs']) * m
    T = np.concatenate([t for t, _ in times])
    Wm = np.zeros((n, len(T)))
    k = 0
    for i, (t, w) in enumerate(times):
        Wm[i, k:k + len(t)] = w
        k += len(t)
    dt = np.abs(T[:, None] - T[None, :])
    rho = Wm @ (p['w'] * np.cos(W12 * dt) + (1 - p['w']) * np.exp(-dt / p['tau'])) @ Wm.T
    K = np.sqrt(np.outer(v, v)) * depth_corr(z, p['l0'], p['l1'], p.get('nu', 1.5), p.get('zl'), p.get('lk')) * rho
    if p.get('u0', 0) > 0:
        B = borehole_vectors(sub)
        cnt = sub['count'].to_numpy(float) if 'count' in sub else np.ones(n)
        t2 = sub['t2'].to_numpy(int)
        Bc = B.copy()
        for pt in np.unique(t2):
            j = t2 == pt
            Bc[j] -= (cnt[j, None] * B[j]).sum(0) / cnt[j].sum()
        u = p['u0'] * np.exp(-z / p.get('zu', np.inf)) if np.isfinite(p.get('zu', np.inf)) else np.full(n, p['u0'])
        K += np.sqrt(np.outer(u, u)) * depth_corr(z, p['lb0'], p['lb1'], p.get('nub', 1.5), p.get('zlb'), p.get('lkb')) * (Bc @ Bc.T)
    if p.get('beta', 0) > 0:                       # slow drift: a misfit trend per entity, smooth in depth
        tm = (Wm @ T) / 12.0                       # row mean time [years]
        tm = tm - tm.mean()
        K += p['beta'] ** 2 * depth_corr(z, p.get('ld', 50.0), 0.0, 1.5) * np.outer(tm, tm)
    if p.get('centre', False):
        cnt = sub['count'].to_numpy(float) if 'count' in sub else np.ones(n)
        t2 = sub['t2'].to_numpy(int)
        Om = np.zeros((n, n))
        for pt in np.unique(t2):
            j = np.flatnonzero(t2 == pt)
            Om[np.ix_(j, j)] = cnt[j] / cnt[j].sum()
        C = np.eye(n) - Om
        K = C @ K @ C.T
    return K


def seasonal_matrix(rows, p, obs=None, times=None):
    times = row_times(rows, obs) if times is None else times
    e = rows['e'].to_numpy(int)
    K = np.zeros((len(rows), len(rows)))
    for ent in np.unique(e):
        idx = np.flatnonzero(e == ent)
        K[np.ix_(idx, idx)] = kernel_block(rows.iloc[idx], [times[i] for i in idx], p)
    return K


def date_spread(rows, obs, Y3_design):
    """Per T3r row: variance of the model over a range row's months (its date is unknown within
    them), from the T3 model values of the same point (one design); 0 for month rows."""
    T3 = obs['T3']
    lut = dict(zip(zip(T3.t2.to_numpy(), T3.ym.to_numpy()), Y3_design))
    out = np.zeros(len(rows))
    for k, (kind, t2, months) in enumerate(zip(rows.kind, rows.t2, rows.months)):
        if kind == 'range':
            out[k] = np.var([lut[(t2, int(mm))] for mm in months])
    return out


def misfit(r, Sigma, rows):
    """r^T Sigma^-1 r / n (entity blocks), the level part (GLS point means, p dims) and contrast
    part (n - p dims), and each row's contribution r_i (Sigma^-1 r)_i."""
    e = rows['e'].to_numpy(int); t2 = rows['t2'].to_numpy(int)
    tot = lvl = 0.0; n_pts = 0
    contrib = np.zeros(len(r))
    per = []
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
        le = u @ np.linalg.solve(A, u)
        lvl += le
        per.append((ent, len(idx), len(pts), rr @ Sir, le))
    n = len(r)
    per = pd.DataFrame(per, columns=['e', 'n', 'p', 'chi2', 'level'])
    return dict(n=n, points=n_pts, ratio=tot / n, chi2=tot, level=lvl / n_pts,
                contrast=(tot - lvl) / max(n - n_pts, 1), contrib=contrib, per_entity=per)
