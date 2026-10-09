"""Stage 2 scratch library: within-point pairs of T3 month rows and their semivariogram."""
import numpy as np, pandas as pd
import target_analysis as ta

LAG_BINS = [(1, 1, '1'), (2, 2, '2'), (3, 3, '3'), (4, 5, '4-5'), (6, 6, '6'), (7, 8, '7-8'),
            (9, 11, '9-11'), (12, 12, '12'), (13, 23, '13-23'), (24, 59, '24-59'),
            (60, 119, '60-119'), (120, 10 ** 6, '>=120')]
BANDS = ['0-5 m', '5-15 m', '15-30 m', '> 30 m']


def lag_bin(lag):
    out = np.full(len(lag), '', dtype=object)
    for lo, hi, name in LAG_BINS:
        out[(lag >= lo) & (lag <= hi)] = name
    return out


def month_rows(obs, mod, common=None):
    """T3 month-dated rows (range_only False) over the common rows: (rows, Y) with Y design x rows."""
    t3 = ta.target('T3', obs, mod, common)
    m = ~t3.rows.range_only.to_numpy()
    rows = t3.rows[m].reset_index(drop=True)
    rows['zb'] = ta.depth_band(rows.depth)
    return rows, t3.Y[:, m], t3.y[m]


def pairs(rows):
    """All pairs (i < j) of month rows of the same pooled point: DataFrame with i, j, t2, e,
    entity, glacier, depth, zb, lag (integer months), dlag (exact months between the readings'
    mean dates), s1, s2 (seasons, sorted), same_bh (one borehole, the same in both rows),
    var_obs (0.5 (sigma_i^2 + sigma_j^2))."""
    out = []
    tm = rows.tmon.to_numpy(int); t = rows.t.to_numpy(float)
    se = rows.season.to_numpy(object); bh = rows.bh_label.to_numpy(object)
    sig = rows.sigma.to_numpy(float)
    for t2, idx in rows.groupby('t2').indices.items():
        if len(idx) < 2:
            continue
        ii, jj = np.triu_indices(len(idx), 1)
        out.append(pd.DataFrame(dict(i=idx[ii], j=idx[jj])))
    P = pd.concat(out, ignore_index=True)
    i, j = P.i.to_numpy(), P.j.to_numpy()
    for c in ('t2', 'e', 'entity', 'glacier', 'system_short', 'depth', 'zb'):
        P[c] = rows[c].to_numpy()[i]
    P['lag'] = np.abs(tm[i] - tm[j])
    P['dlag'] = np.abs(t[i] - t[j]) * 12.0
    a, b = se[i], se[j]
    order = np.array([ta.SEASONS.index(x) for x in a]) <= np.array([ta.SEASONS.index(x) for x in b])
    P['s1'] = np.where(order, a, b); P['s2'] = np.where(order, b, a)
    P['same_bh'] = (bh[i] == bh[j]) & ~pd.Series(bh[i]).str.contains('/').to_numpy()
    P['var_obs'] = 0.5 * (sig[i] ** 2 + sig[j] ** 2)
    P['lagbin'] = lag_bin(P.lag.to_numpy())
    return P


def gamma(P, R):
    """Half squared differences of residuals R (n_sets x rows) for pairs P: (n_sets x pairs)."""
    R = np.atleast_2d(R)
    return 0.5 * (R[:, P.i.to_numpy()] - R[:, P.j.to_numpy()]) ** 2


def binned(P, G, by, balanced=False):
    """Mean gamma per group of `by` columns for each set (rows of G): DataFrame (groups x sets).
    balanced: each pooled point weighs equally within a group (mean of its pairs first)."""
    keys = [P[c].to_numpy() for c in by]
    df = pd.DataFrame(G.T)
    for c, k in zip(by, keys):
        df[c] = k
    if balanced:
        df['t2'] = P.t2.to_numpy()
        df = df.groupby(by + ['t2']).mean().groupby(level=list(range(len(by)))).mean()
    else:
        df = df.groupby(by).mean()
    return df
