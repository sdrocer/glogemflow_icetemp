"""Misfit ratios r^T Sigma^-1 r / n with Sigma = D + S + K (+ date spread for T3r range rows)."""
import sys, pickle, time; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2kern2 as K2
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
pw = pickle.load(open(W + 'within_d13.pkl', 'rb'))['exp'][0]
rest = pickle.load(open(W + 'rest_d13.pkl', 'rb'))
PA = dict(a=pw['a'], zs=pw['zs'], w=pw['w'], tau=pw['tau'], q=pw['q'], nu=1.5,
          l0=rest[('within', 1.5)][0], l1=rest[('within', 1.5)][1],
          u0=rest[('u', 'const')][0], lb0=0.05, lb1=rest[('between', 1.5)][1], nub=1.5)
if len(sys.argv) > 1:
    PA.update(eval(sys.argv[1]))
print('kernel:', {k: round(float(v), 4) for k, v in PA.items()})
# spike rows (T3 month rows) -> keys (t2, ym)
mrows = ta.target('T3', obs, mod).rows
mrows = mrows[~mrows.range_only].reset_index(drop=True)
spike = np.load(W + 'spike_rows_month.npy')
spike_keys = set(zip(mrows.t2.to_numpy()[spike], mrows.ym.to_numpy()[spike]))
y2 = ta.target('T2', obs, mod)
rm2 = np.sqrt(np.mean((y2.Y - y2.y) ** 2, axis=1))
best5 = list(np.argsort(rm2)[:5])
med_run = int(np.argsort(rm2)[len(rm2) // 2])
print('best five by T2 RMSE:', best5, 'median run:', med_run, f'(T2 RMSE {rm2[med_run]:.3f})')
runs = [13] + [b for b in best5 if b != 13] + [med_run]
out = []
for name in ('T2', 'T3', 'T3r', 'T3m'):
    t = ta.target('T3' if name == 'T3m' else name, obs, mod)
    rows, Y, y = t.rows, t.Y, t.y
    if name == 'T3m':
        m = ~rows.range_only.to_numpy()
        rows, Y, y = rows[m].reset_index(drop=True), Y[:, m], y[m]
    t0 = time.time()
    Kmat = K2.seasonal_matrix(rows, PA, obs=obs)
    base = ta.error_cov(rows)
    ym_key = rows.ym.to_numpy() if 'ym' in rows else (rows.ym_first.to_numpy() if 'kind' in rows else None)
    is_spike = np.array([(a, b) in spike_keys for a, b in zip(rows.t2, ym_key)]) if ym_key is not None else np.zeros(len(rows), bool)
    if 'kind' in rows:
        is_spike &= (rows.kind == 'month').to_numpy()
    for run in runs:
        ds = K2.date_spread(rows, obs, mod['Y3'][run]) if name == 'T3r' else np.zeros(len(rows))
        r = Y[run] - y
        for lab, Sig in (('D+S', base), ('D+S+K', base + Kmat + np.diag(ds))):
            mm = K2.misfit(r, Sig, rows)
            rec = dict(target=name, run=run, sigma=lab, n=mm['n'], ratio=mm['ratio'], level=mm['level'], contrast=mm['contrast'])
            if lab == 'D+S+K':
                keep = ~is_spike
                if is_spike.any():
                    mm2 = K2.misfit(r[keep], Sig[np.ix_(keep, keep)], rows[keep].reset_index(drop=True))
                    rec.update(ratio_nospike=mm2['ratio'], contrast_nospike=mm2['contrast'], n_spike=int(is_spike.sum()))
                c = pd.Series(mm['contrib']).groupby(rows.system_short.to_numpy()).agg(['sum', 'size'])
                for s_, row_ in c.iterrows():
                    rec[f'sys {s_}'] = row_['sum'] / row_['size']
                if run == 13:
                    pickle.dump(dict(rows=rows, contrib=mm['contrib'], r=r), open(W + f'contrib_{name}_d13.pkl', 'wb'))
            out.append(rec)
    print(f'{name}: {len(rows)} rows ({time.time()-t0:.0f} s)', flush=True)
O = pd.DataFrame(out)
pd.set_option('display.float_format', lambda v: f'{v:.3f}')
print(O[O.run == 13].drop(columns='run').to_string(index=False))
print('\nall runs evaluated (ratio with D+S+K):')
print(O[O.sigma == 'D+S+K'].pivot(index='run', columns='target', values='ratio').loc[runs].to_string())
O.to_pickle(W + 'ratios_A.pkl')
