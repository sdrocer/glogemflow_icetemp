"""Part 1: within-borehole month-to-month semivariogram (same-borehole pairs, spike rows out)."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
from scipy.optimize import minimize
import target_analysis as ta, s2lib as L
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
W12 = 2 * np.pi / 12
KN = np.array([1.0, 4.0, 10.0, 20.0, 40.0])
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(W + 'spike_rows_month.npy')
P0 = pd.read_pickle(W + 'pairs.pkl')
P = P0[~(P0.i.isin(spike) | P0.j.isin(spike)) & P0.same_bh].reset_index(drop=True)
ok = np.ones(len(rows), bool); ok[spike] = False
# contrast weights within each (point, borehole) series: 2 / rows of that series
rows['bh1'] = rows.bh_label
ser = rows[ok].groupby(['t2', 'bh1']).size()
P['nser'] = [ser[(t, b)] for t, b in zip(P.t2, rows.bh1.to_numpy()[P.i])]
P['w'] = 2.0 / P.nser
P['ice'] = rows.firn_label.to_numpy()[P.i] == 'ice'
R = Y - y[None, :]


def vz(p, z, form):
    if form == 'exp':
        return p['a'] * np.exp(-z / p['zs'])
    return np.exp(np.interp(z, KN, [p[f'v{k}'] for k in range(5)]))


def gmod(p, z, lag, ice, form):
    v = vz(p, z, form) * np.where(ice, p['q'], 1.0)
    rho = p['w'] * np.cos(W12 * lag) + (1 - p['w']) * np.exp(-lag / p['tau'])
    return v * (1 - rho)


def fit(r, form, start=None):
    d2 = (r[P.i] - r[P.j]) ** 2
    z, lag, ice, vobs, w = P.depth.to_numpy(), P.dlag.to_numpy(), P.ice.to_numpy(), P.var_obs.to_numpy(), P.w.to_numpy()
    if form == 'exp':
        names = ['a', 'zs', 'w', 'tau', 'q']; lo = [-9, np.log(0.2), 0, np.log(0.05), np.log(1e-3)]; hi = [np.log(50), np.log(300), 1, np.log(600), np.log(5)]
        x0 = [np.log(1.5), np.log(5), 0.4, np.log(2), np.log(0.3)]; isl = [1, 1, 0, 1, 1]
    else:
        names = [f'v{k}' for k in range(5)] + ['w', 'tau', 'q']
        lo = [np.log(1e-5)] * 5 + [0, np.log(0.05), np.log(1e-3)]; hi = [np.log(50)] * 5 + [1, np.log(600), np.log(5)]
        x0 = list(np.log([1.0, 0.8, 0.3, 0.05, 0.005])) + [0.4, np.log(2), np.log(0.3)]; isl = [0] * 5 + [0, 1, 1]
    isl = np.array(isl, bool)
    if start is not None:
        x0 = start

    def unpack(x):
        x = np.clip(x, lo, hi)
        return dict(zip(names, np.where(isl, np.exp(x), x)))

    def f(x):
        g = gmod(unpack(x), z, lag, ice, form) + vobs
        return np.sum(w * (0.5 * np.log(2 * g) + d2 / (4 * g)))
    best = None
    for _ in range(4):
        res = minimize(f, x0 if best is None else best.x, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                       options=dict(maxiter=4000, ftol=1e-14, gtol=1e-10))
        res2 = minimize(f, res.x, method='Nelder-Mead', options=dict(maxfev=8000, xatol=1e-7, fatol=1e-9, adaptive=True))
        c = res2 if res2.fun < res.fun else res
        if best is None or c.fun < best.fun - 1e-9:
            best = c
        else:
            break
    return unpack(best.x), -best.fun, best.x


if __name__ == '__main__':
    out = {}
    for form in ('exp', 'knots'):
        p, cl, x = fit(R[13], form)
        out[form] = (p, cl, x)
        print(form, f'CL {cl:.2f}', {k: round(float(v), 4) for k, v in p.items()})
    pickle.dump(out, open(W + 'within_d13.pkl', 'wb'))
    lb = [(1, 1, '1'), (2, 2, '2'), (3, 3, '3'), (4, 5, '4-5'), (6, 6, '6'), (7, 8, '7-8'), (9, 11, '9-11'), (12, 12, '12'),
          (13, 23, '13-23'), (24, 36, '24-36'), (37, 119, '37-119')]
    P['lb'] = ''
    for lo_, hi_, n in lb:
        P.loc[(P.lag >= lo_) & (P.lag <= hi_), 'lb'] = n
    names = [n for *_, n in lb]
    P['meas'] = 0.5 * (R[13, P.i] - R[13, P.j]) ** 2
    P['type'] = np.where(P.ice, 'ice', 'firn')
    for form in ('exp', 'knots'):
        P['mod'] = gmod(out[form][0], P.depth.to_numpy(), P.dlag.to_numpy(), P.ice.to_numpy(), form) + P.var_obs
        g = P.groupby(['type', 'zb', 'lb'])
        t = pd.DataFrame({'meas': g.apply(lambda s: np.sum(s.w * s.meas) / s.w.sum()),
                          'model': g.apply(lambda s: np.sum(s.w * s['mod']) / s.w.sum()), 'pts': g.t2.nunique()})
        print(f'\n== {form}: contrast-weighted measured vs fitted (same borehole) ==')
        print(t.unstack('lb').reindex(columns=names, level=1).round(3).to_string())
        P['ratio'] = P.meas / P['mod']
        print('measured/model, contrast-weighted:', (P.assign(wr=P.w * P.ratio).groupby(['type', 'zb']).apply(lambda s: s.wr.sum() / s.w.sum())).round(2).to_dict())
        print('by system:', (P.assign(wr=P.w * P.ratio).groupby('system_short').apply(lambda s: s.wr.sum() / s.w.sum())).round(2).to_dict())
