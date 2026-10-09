import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2lib as L, s2vfit3 as V
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
spike = np.load(W + 'spike_rows_month.npy')
P0 = pd.read_pickle(W + 'pairs.pkl')
P = P0[~(P0.i.isin(spike) | P0.j.isin(spike))].reset_index(drop=True)
ok = np.ones(len(rows), bool); ok[spike] = False
cnt = rows[ok].groupby('t2').size()
w = (2.0 / P.t2.map(cnt)).to_numpy()
ice = (rows.firn_label.to_numpy()[P.i] == 'ice')
z, lag, vobs = P.depth.to_numpy(), P.dlag.to_numpy(), P.var_obs.to_numpy()
r = Y[13] - y
d2 = (r[P.i] - r[P.j]) ** 2
data = (z, lag, ice, d2, vobs, w)
lv = dict(v0=np.log(1.0), v1=np.log(0.8), v2=np.log(0.3), v3=np.log(0.05), v4=np.log(0.005))
C = {
 'K1 knots, w, tau':                 dict(lv, w=0.4, tau=3.0),
 'K2 + ice factor':                  dict(lv, w=0.4, tau=3.0, q=0.3),
 'K3 + drift':                       dict(lv, w=0.4, tau=3.0, q=0.3, beta=0.05),
 'K4 + drift ice factor':            dict(lv, w=0.4, tau=3.0, q=0.3, beta=0.05, qd=0.3),
 'K5 K3 + interannual':              dict(lv, w=0.4, tau=3.0, q=0.3, beta=0.05, ay=0.5, ty=36.0),
}
fits = {}
for name, st in C.items():
    p, cl, k = V.fit(st, data)
    fits[name] = (p, cl, k)
    vk = ' '.join(f'{np.exp(p[f"v{j}"]):.3g}' for j in range(5))
    oth = ' '.join(f'{a}={b:.4g}' for a, b in p.items() if not a.startswith('v'))
    print(f'{name:28s} CL {cl:9.2f} k {k:2d} AIC {-2*cl+2*k:8.2f}  v(1,4,10,20,40 m)= {vk}  {oth}')
pickle.dump(dict(fits=fits), open(W + 'vfit3_d13.pkl', 'wb'))
P['meas'] = d2 / 2; P['w'] = w; P['type'] = np.where(ice, 'ice', 'firn')
lb = [(1, 2, '1-2'), (3, 5, '3-5'), (6, 8, '6-8'), (9, 11, '9-11'), (12, 12, '12'), (13, 36, '13-36'), (37, 119, '37-119'), (120, 10**6, '>=120')]
P['lb'] = ''
for lo, hi, n in lb:
    P.loc[(P.lag >= lo) & (P.lag <= hi), 'lb'] = n
names = [n for *_, n in lb]
for name in ('K3 + drift', 'K5 K3 + interannual'):
    P['mod'] = V.gamma_model(fits[name][0], z, lag, ice) + vobs
    g = P.groupby(['type', 'zb', 'lb'])
    t = pd.DataFrame({'meas': g.apply(lambda s: np.sum(s.w * s.meas) / s.w.sum()),
                      'model': g.apply(lambda s: np.sum(s.w * s['mod']) / s.w.sum())})
    print(f'\n== {name}: contrast-weighted measured vs fitted ==')
    print(t.unstack('lb').reindex(columns=names, level=1).round(3).to_string())
    # aggregate ratio of measured to model, contrast-weighted, by type x band
    P['ratio'] = P.meas / P['mod']
    print((P.assign(wr=P.w * P.ratio).groupby(['type', 'zb']).apply(lambda s: s.wr.sum() / s.w.sum())).round(2).to_string())
