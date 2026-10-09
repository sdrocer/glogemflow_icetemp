"""Depth correlation of the within-borehole part (same-borehole increments at two depths over the
same two months) and of the between-borehole part (differences between the same two boreholes at
two depths): measured binned correlations, and the length l(z) = l0 + l1 z of a non-stationary
exponential correlation fitted to them (least squares over bins, weight = increment pairs)."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(W + 'spike_rows_month.npy')
R = Y - y[None, :]
P0 = pd.read_pickle(W + 'pairs.pkl')
P = P0[~(P0.i.isin(spike) | P0.j.isin(spike))].reset_index(drop=True)
P['bh_i'] = rows.bh_label.to_numpy()[P.i]; P['bh_j'] = rows.bh_label.to_numpy()[P.j]
P['ym_i'] = rows.ym.to_numpy()[P.i]; P['ym_j'] = rows.ym.to_numpy()[P.j]
P['ice'] = rows.firn_label.to_numpy()[P.i] == 'ice'
pw = pickle.load(open(W + 'within_d13.pkl', 'rb'))['exp'][0]


def lengths(z, l0, l1):
    return l0 + l1 * z


def pac(z1, z2, l0, l1):
    a, b = lengths(z1, l0, l1), lengths(z2, l0, l1)
    L2 = 0.5 * (a * a + b * b)
    return np.sqrt(a * b / L2) * np.exp(-np.abs(z1 - z2) / np.sqrt(L2))


def build(same):
    Q = P[P.same_bh == same]
    recs = []
    for k, g in Q.groupby(['e', 'bh_i', 'bh_j', 'ym_i', 'ym_j']):
        if len(g) < 2:
            continue
        a, b = np.triu_indices(len(g), 1)
        zz, ii, jj = g.depth.to_numpy(), g.i.to_numpy(), g.j.to_numpy()
        lag, ice = g.dlag.to_numpy(), g.ice.to_numpy()
        for x, u in zip(a, b):
            recs.append((k[0], zz[x], zz[u], ii[x], jj[x], ii[u], jj[u], lag[x], ice[x], ice[u]))
    C = pd.DataFrame(recs, columns=['e', 'z1', 'z2', 'i1', 'j1', 'i2', 'j2', 'lag', 'ice1', 'ice2'])
    C['zmin'] = np.minimum(C.z1, C.z2); C['dz'] = np.abs(C.z1 - C.z2)
    C['zb'] = pd.cut(C.zmin, [0, 5, 15, 30, 200], labels=['0-5', '5-15', '15-30', '>30'])
    C['dzb'] = pd.cut(C.dz, [0, 1.5, 2.5, 3.5, 5, 7.5, 10, 15, 25, 200])
    return C


def measured(C, run=13):
    D1 = R[run, C.i1] - R[run, C.j1]; D2 = R[run, C.i2] - R[run, C.j2]
    t = pd.DataFrame(dict(p=D1 * D2, s1=D1 ** 2, s2=D2 ** 2, zb=C.zb, dzb=C.dzb))
    g = t.groupby(['zb', 'dzb'], observed=True)
    return pd.DataFrame(dict(n=g.size(), corr=g.p.sum() / np.sqrt(g.s1.sum() * g.s2.sum())))


def model_within(C, l0, l1):
    v1 = pw['a'] * np.exp(-C.z1 / pw['zs']) * np.where(C.ice1, pw['q'], 1)
    v2 = pw['a'] * np.exp(-C.z2 / pw['zs']) * np.where(C.ice2, pw['q'], 1)
    one_m_rho = 1 - (pw['w'] * np.cos(W12 * C.lag) + (1 - pw['w']) * np.exp(-C.lag / pw['tau']))
    cov = 2 * np.sqrt(v1 * v2) * pac(C.z1.to_numpy(), C.z2.to_numpy(), l0, l1) * one_m_rho
    return cov, 2 * v1 * one_m_rho, 2 * v2 * one_m_rho


def model_between(C, l0, l1):   # constant variance: correlation = depth correlation
    c = pac(C.z1.to_numpy(), C.z2.to_numpy(), l0, l1)
    return c, np.ones(len(C)), np.ones(len(C))


def binned_model(C, fn, l0, l1):
    cov, v1, v2 = fn(C, l0, l1)
    t = pd.DataFrame(dict(c=cov, v1=v1, v2=v2, zb=C.zb, dzb=C.dzb))
    g = t.groupby(['zb', 'dzb'], observed=True)
    return g.c.sum() / np.sqrt(g.v1.sum() * g.v2.sum())


def fit(C, fn, M):
    def f(x):
        l0, l1 = np.exp(x[0]), x[1]
        if l1 < 0:
            return 1e9
        mc = binned_model(C, fn, l0, l1)
        return float(np.sum(M['n'] * (M['corr'] - mc) ** 2))
    best = min((minimize(f, [np.log(s0), s1], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxfev=4000))
                for s0 in (0.5, 2, 8) for s1 in (0.0, 0.3, 1.0)), key=lambda r: r.fun)
    return np.exp(best.x[0]), best.x[1], best.fun


out = {}
for label, same, fn in (('within', True, model_within), ('between', False, model_between)):
    C = build(same)
    M = measured(C)
    Mmed = pd.concat([measured(C, r)['corr'] for r in range(32)], axis=1).median(axis=1)
    l0, l1, sse = fit(C, fn, M)
    lc0, _, sse_c = fit(C, lambda C_, a, b: fn(C_, a, 0.0), M)
    M['corr_med'] = Mmed
    M['model'] = binned_model(C, fn, l0, l1)
    M['model_const'] = binned_model(C, fn, lc0, 0.0)
    out[label] = dict(l0=l0, l1=l1, sse=sse, l_const=lc0, sse_const=sse_c, table=M)
    print(f'\n== {label}-borehole increments: {len(C)} pairs; fitted l(z) = {l0:.3f} + {l1:.3f} z (SSE {sse:.1f}); '
          f'constant l = {lc0:.2f} m (SSE {sse_c:.1f}) ==')
    print(M.round(3).to_string())
pickle.dump(out, open(W + 'depthcorr_d13.pkl', 'wb'))
