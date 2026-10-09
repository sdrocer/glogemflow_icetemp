"""Stage 2 report numbers, from target_analysis (module) and the Stage 2 scratch pairs.
Writes stage2/report.log and stage2/report.json."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import sys, json, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2lib as L
pd.set_option('display.width', 250, 'display.max_columns', 80, 'display.max_rows', 500)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
p = ta.SEASONAL
keep29 = mod['design'].keep.to_numpy()
rows, Y, y = L.month_rows(obs, mod)                 # T3 month-dated rows over the common rows
R = Y - y[None, :]
spike = np.load(W + 'spike_rows_month.npy')
P0 = pd.read_pickle(W + 'pairs.pkl')                # all same-point pairs of month rows
P0['spike'] = (P0.i.isin(spike) | P0.j.isin(spike)).to_numpy()
ok = np.ones(len(rows), bool); ok[spike] = False
npt = rows[ok].groupby('t2').size()
P = P0[~P0.spike].reset_index(drop=True)
P['w'] = 2.0 / P.t2.map(npt).to_numpy()            # each point's independent differences count once
bh = rows.bh_label.to_numpy()


def bdist(a, b):
    a, b = set(str(a).split('/')), set(str(b).split('/'))
    k = sorted(a | b)
    x = np.array([1 / len(a) if q in a else 0 for q in k]); yv = np.array([1 / len(b) if q in b else 0 for q in k])
    return float(np.sum((x - yv) ** 2))


P['bd'] = [bdist(bh[i], bh[j]) for i, j in zip(P.i, P.j)]
P['ice'] = rows.firn_label.to_numpy()[P.i] == 'ice'
P['high'] = (rows.band.to_numpy()[P.i] >= p['high_band']) & ~P.ice
P['kind'] = np.where(P.ice, 'ice', np.where(P.high, 'firn >= 4000 m', 'firn < 4000 m'))
m = np.where(P.ice, p['ice'], np.where(P.high, p['high'], 1.0))
lag = P.dlag.to_numpy()
P['fit'] = (m * p['a'] * np.exp(-P.depth / p['zs']) * (p['w'] * (1 - np.cos(2 * np.pi * lag / 12)) + (1 - p['w']) * (1 - np.exp(-lag / p['tau'])))
            + p['u0'] * P.bd / 2 + 0.5 * (p['beta'] * lag / 12) ** 2 + P.var_obs)
G = 0.5 * (R[:, P.i] - R[:, P.j]) ** 2              # runs x pairs
lagnames = [n for *_, n in L.LAG_BINS]
log = []


def emit(s=''):
    print(s); log.append(str(s))


def table(Q, Gq, by, label):
    """Contrast-weighted mean gamma per group: design 13, median over runs (32 and 29), fitted."""
    keys = [Q[c].to_numpy() for c in by]
    idx = pd.MultiIndex.from_arrays(keys, names=by)
    wq = Q.w.to_numpy()
    num = pd.DataFrame(Gq.T * wq[:, None], index=idx).groupby(level=list(range(len(by)))).sum()
    den = pd.Series(wq, index=idx).groupby(level=list(range(len(by)))).sum()
    per_run = num.div(den, axis=0)
    out = pd.DataFrame({'d13': per_run[13], 'med32': per_run.median(axis=1), 'med29': per_run.loc[:, np.flatnonzero(keep29)].median(axis=1),
                        'fit': pd.Series(Q.fit.to_numpy() * wq, index=idx).groupby(level=list(range(len(by)))).sum() / den,
                        'pairs': pd.Series(1, index=idx).groupby(level=list(range(len(by)))).sum(),
                        'points': pd.Series(Q.t2.to_numpy(), index=idx).groupby(level=list(range(len(by)))).nunique()})
    emit(f'\n== {label} =='); return out


res = {}
# 1. semivariogram, within one borehole (the seasonal part)
S1 = P[P.bd == 0]
t = table(S1, G[:, (P.bd == 0).to_numpy()], ['zb', 'lagbin'], 'semivariogram within one borehole [C^2], contrast-weighted (spike rows out)')
tw = t.unstack('zb').reindex(lagnames)
emit(tw[['d13', 'med32', 'fit', 'points']].reindex(columns=L.BANDS, level=1).round(3).to_string())
res['within_borehole'] = json.loads(tw.round(4).to_json())
# 1b. by kind (firn < 4000 m, firn >= 4000 m, ice) at grouped lags
lb = [(1, 2, '1-2'), (3, 5, '3-5'), (6, 8, '6-8'), (9, 11, '9-11'), (12, 12, '12'), (13, 36, '13-36'), (37, 10 ** 6, '>36')]
P['lb'] = ''
for lo, hi, n in lb:
    P.loc[(P.lag >= lo) & (P.lag <= hi), 'lb'] = n
S1 = P[P.bd == 0]
t = table(S1, G[:, (P.bd == 0).to_numpy()], ['kind', 'zb', 'lb'], 'within one borehole by site kind (contrast-weighted)')
tk = t[['d13', 'med32', 'fit', 'points']].unstack('lb').reindex(columns=[n for *_, n in lb], level=1)
emit(tk.round(3).to_string())
res['within_by_kind'] = json.loads(t.round(4).reset_index().to_json(orient='records'))
# 2. between boreholes of one point
S2 = P[P.bd > 0].copy()
S2['yb'] = pd.cut(S2.dlag / 12, [0, 5, 10, 20, 50], labels=['<5 yr', '5-10 yr', '10-20 yr', '20-50 yr'])
t = table(S2, G[:, (P.bd > 0).to_numpy()], ['zb', 'yb'], 'between boreholes of one point (all Monte Rosa) [C^2]')
emit(t.unstack('yb')[['d13', 'med32', 'fit', 'points']].reindex(L.BANDS).round(3).to_string())
res['between_boreholes'] = json.loads(t.round(4).reset_index().astype({'yb': str}).to_json(orient='records'))
# 3. season pairs (within one borehole, lag <= 12 months)
S3 = P[(P.bd == 0) & (P.lag <= 12)]
t = table(S3, G[:, ((P.bd == 0) & (P.lag <= 12)).to_numpy()], ['zb', 's1', 's2'], 'season pairs, one borehole, lag <= 12 months')
t['mean_lag'] = S3.groupby(['zb', 's1', 's2']).lag.mean()
emit(t.loc[['0-5 m', '5-15 m']].round(3).to_string())
res['season_pairs'] = json.loads(t.round(4).reset_index().to_json(orient='records'))
# 4. level vs month-to-month (points with >= 2 month rows, spike rows out)
Rk = R[:, ok]; rk = rows[ok].reset_index(drop=True)
grp = rk.groupby('t2').indices
multi = {k: v for k, v in grp.items() if len(v) >= 2}
recs = []
for k, ix in multi.items():
    rr = Rk[:, ix]
    recs.append(dict(t2=k, zb=rk.zb.iloc[ix[0]], n=len(ix), level=rr.mean(axis=1) ** 2,
                     m2m=rr.var(axis=1, ddof=1)))
Dl = pd.DataFrame(recs)
lv = np.vstack(Dl.level.to_numpy()); mv = np.vstack(Dl.m2m.to_numpy())
split = pd.DataFrame({'points': Dl.groupby('zb').size(),
                      'level d13': pd.Series(lv[:, 13]).groupby(Dl.zb.to_numpy()).mean(),
                      'level med32': pd.Series(np.median(lv, axis=1)).groupby(Dl.zb.to_numpy()).mean(),
                      'm2m d13': pd.Series(mv[:, 13]).groupby(Dl.zb.to_numpy()).mean(),
                      'm2m med32': pd.Series(np.median(mv, axis=1)).groupby(Dl.zb.to_numpy()).mean()}).reindex(L.BANDS)
emit('\n== level (point-mean misfit squared) vs month-to-month (within-point variance) [C^2], points with >= 2 month rows ==')
emit(split.round(3).to_string())
res['level_vs_m2m'] = json.loads(split.round(4).to_json())
# 5. misfit ratios
y2 = ta.target('T2', obs, mod)
rm2 = np.sqrt(np.mean((y2.Y - y2.y) ** 2, axis=1))
order = np.argsort(rm2)
best5, med_run = [int(x) for x in order[:5]], int(order[len(order) // 2])
mk = {}
spike_keys = set(zip(rows.t2.to_numpy()[spike], rows.ym.to_numpy()[spike]))


def tgt_rows(name):
    t = ta.target('T3' if name == 'T3 month rows' else name, obs, mod)
    rr, Yt, yt = t.rows, t.Y, t.y
    if name == 'T3 month rows':
        k = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[k].reset_index(drop=True), Yt[:, k], yt[k]
    return t, rr, Yt, yt


rat = []
for name in ('T1', 'T2', 'T3', 'T3r', 'T3 month rows'):
    t, rr, Yt, yt = tgt_rows(name)
    if name == 'T1':
        C = ta.target_error_cov(t, obs)
        for run in [13] + best5[1:] + [med_run]:
            r = Yt[run] - yt
            rat.append(dict(target=name, error='D+S+K', run=run, ratio=float(r @ np.linalg.solve(C, r)) / len(r), level=np.nan, m2m=np.nan, n=len(r)))
        continue
    base = ta.error_cov(rr)
    K = ta.seasonal_cov(rr)
    mm = ta.misfit_ratio(Yt[13] - yt, base, rr)
    rat.append(dict(target=name, error='D+S', run=13, ratio=mm['ratio'], level=mm['level'], m2m=mm['month_to_month'], n=mm['n']))
    key = rr.ym.to_numpy() if 'ym' in rr else (rr.ym_first.to_numpy() if 'kind' in rr else None)
    sp = np.array([(a, b) in spike_keys for a, b in zip(rr.t2, key)]) if key is not None else np.zeros(len(rr), bool)
    if 'kind' in rr:
        sp &= (rr.kind == 'month').to_numpy()
    for run in [13] + best5[1:] + [med_run]:
        C = base + K + (np.diag(ta.date_spread_var(rr, obs, mod['Y3'][run])) if name == 'T3r' else 0.0)
        r = Yt[run] - yt
        mm = ta.misfit_ratio(r, C, rr)
        rec = dict(target=name, error='D+S+K', run=run, ratio=mm['ratio'], level=mm['level'], m2m=mm['month_to_month'], n=mm['n'])
        if sp.any():
            k = ~sp
            m2 = ta.misfit_ratio(r[k], C[np.ix_(k, k)], rr[k].reset_index(drop=True))
            rec.update(ratio_no_spike=m2['ratio'], m2m_no_spike=m2['month_to_month'])
        rat.append(rec)
RT = pd.DataFrame(rat)
emit(f'\n== misfit ratio r^T Sigma^-1 r / n (best five runs by T2 RMSE {best5}, median run {med_run}) ==')
emit(RT.round(3).to_string(index=False))
res['ratios'] = json.loads(RT.round(4).to_json(orient='records'))
res['runs'] = dict(best5=best5, median_run=med_run, t2_rmse=dict(zip(map(str, range(32)), np.round(rm2, 4).tolist())))
# 6. T3r at design 13 by depth band and by system: level and month-to-month parts
t, rr, Yt, yt = tgt_rows('T3r')
C = ta.error_cov(rr) + ta.seasonal_cov(rr) + np.diag(ta.date_spread_var(rr, obs, mod['Y3'][13]))
r = Yt[13] - yt
by = []
for lab, mask in ([('depth ' + b, ta.depth_band(rr.depth) == b) for b in L.BANDS]
                  + [('system ' + s, (rr.system_short == s).to_numpy()) for s in sorted(rr.system_short.unique())]):
    # restrict to whole points: a point's rows share depth and system, so masks keep points whole
    sub = rr[mask].reset_index(drop=True)
    mm = ta.misfit_ratio(r[mask], C[np.ix_(mask, mask)], sub)
    by.append(dict(group=lab, rows=mm['n'], points=mm['points'], ratio=mm['ratio'], level=mm['level'], m2m=mm['month_to_month']))
BY = pd.DataFrame(by)
emit('\n== T3r design 13 by depth band and system (each group on its own; cross-group correlations within an entity dropped) ==')
emit(BY.round(3).to_string(index=False))
res['t3r_by_group'] = json.loads(BY.round(4).to_json(orient='records'))
# 7. sensitivity of T3r (design 13)
sens = []
for lab, ch in [('proposed', {}), ('no high-firn factor', dict(high=1.0)), ('no between-borehole part', dict(u0=0.0)),
                ('no drift', dict(beta=0.0)), ('periodic share 0', dict(w=0.0)), ('periodic share 0.5', dict(w=0.5)),
                ('depth lengths x0.5', dict(lens=tuple(0.5 * np.array(p['lens'])))), ('depth lengths x2', dict(lens=tuple(2 * np.array(p['lens'])))),
                ('variance a x0.5', dict(a=p['a'] * 0.5)), ('variance a x2', dict(a=p['a'] * 2))]:
    Kc = ta.seasonal_cov(rr, ch)
    mm = ta.misfit_ratio(r, ta.error_cov(rr) + Kc + np.diag(ta.date_spread_var(rr, obs, mod['Y3'][13])), rr)
    sens.append(dict(variant=lab, ratio=mm['ratio'], level=mm['level'], m2m=mm['month_to_month']))
SE = pd.DataFrame(sens)
emit('\n== T3r sensitivity (design 13) =='); emit(SE.round(3).to_string(index=False))
res['sensitivity'] = json.loads(SE.round(4).to_json(orient='records'))
# 8. fitted parameters per run
T = pd.read_pickle(W + 'final_fit_runs.pkl')
pr = pd.DataFrame({'design 13': T.loc[13], 'median 29': T[keep29].median(), 'min 29': T[keep29].min(), 'max 29': T[keep29].max()}).drop('cl')
emit('\n== fitted parameters (pair likelihood), design 13 and over the 29 runs with f_ref >= 0.1 =='); emit(pr.round(4).to_string())
res['fit_per_run'] = json.loads(pr.round(4).to_json())
res['proposed'] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in p.items()}
res['spike_rows'] = rows.iloc[spike][['entity', 'depth', 'ym', 'bh_label', 'obs']].assign(model_d13=Y[13, spike]).round(2).to_dict(orient='records')
open(W + 'report.log', 'w').write('\n'.join(log) + '\n')
json.dump(res, open(W + 'report.json', 'w'), indent=1, default=float)
