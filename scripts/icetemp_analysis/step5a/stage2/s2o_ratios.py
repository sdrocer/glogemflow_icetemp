"""Misfit ratios with the final kernel (Sigma = D + S + K, plus date spread for T3r range rows)."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta, s2kern2 as K2
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 400)
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
FINAL = pickle.load(open(W + 'final_kernel.pkl', 'rb'))['final']
mrows = ta.target('T3', obs, mod).rows
mrows = mrows[~mrows.range_only].reset_index(drop=True)
spike = np.load(W + 'spike_rows_month.npy')
spike_keys = set(zip(mrows.t2.to_numpy()[spike], mrows.ym.to_numpy()[spike]))
y2 = ta.target('T2', obs, mod)
rm2 = np.sqrt(np.mean((y2.Y - y2.y) ** 2, axis=1))
order = np.argsort(rm2)
runs = [int(x) for x in order[:5]] + [int(order[len(order) // 2])]
TG = {}
for name in ('T2', 'T3', 'T3r', 'T3m'):
    t = ta.target('T3' if name == 'T3m' else name, obs, mod)
    rr, Yt, yt = t.rows, t.Y, t.y
    if name == 'T3m':
        m = ~rr.range_only.to_numpy(); rr, Yt, yt = rr[m].reset_index(drop=True), Yt[:, m], yt[m]
    key = rr.ym.to_numpy() if 'ym' in rr else (rr.ym_first.to_numpy() if 'kind' in rr else None)
    sp = np.array([(a, b) in spike_keys for a, b in zip(rr.t2, key)]) if key is not None else np.zeros(len(rr), bool)
    if 'kind' in rr:
        sp &= (rr.kind == 'month').to_numpy()
    TG[name] = (rr, Yt, yt, sp)


def evaluate(p, run, names=('T2', 'T3', 'T3r', 'T3m'), by_system=False):
    out = {}
    for name in names:
        rr, Yt, yt, sp = TG[name]
        Sig = ta.error_cov(rr) + (K2.seasonal_matrix(rr, p, obs=obs) if p is not None else 0.0)
        if name == 'T3r' and p is not None:
            Sig = Sig + np.diag(K2.date_spread(rr, obs, mod['Y3'][run]))
        r = Yt[run] - yt
        mm = K2.misfit(r, Sig, rr)
        rec = dict(ratio=mm['ratio'], level=mm['level'], contrast=mm['contrast'], n=mm['n'], points=mm['points'])
        if sp.any():
            k = ~sp
            m2 = K2.misfit(r[k], Sig[np.ix_(k, k)], rr[k].reset_index(drop=True))
            rec.update(ratio_nospike=m2['ratio'], contrast_nospike=m2['contrast'])
        if by_system:
            pe = mm['per_entity']
            pe['system'] = rr.groupby('e').system_short.first().reindex(pe.e).to_numpy()
            g = pe.groupby('system').agg(n=('n', 'sum'), p=('p', 'sum'), chi2=('chi2', 'sum'), level=('level', 'sum'))
            g['level_ratio'] = g.level / g.p
            g['contrast_dims'] = g.n - g.p
            g['contrast_ratio'] = (g.chi2 - g.level) / g.contrast_dims.clip(lower=1)
            g['ratio'] = g.chi2 / g.n
            rec['by_system'] = g
            rec['contrib'] = mm['contrib']
        out[name] = rec
    return out


if __name__ == '__main__':
    res = {}
    for run in runs:
        res[run] = evaluate(FINAL, run, by_system=(run == 13))
    base13 = evaluate(None, 13)
    print('kernel', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in FINAL.items()})
    print(f'runs: best five by T2 RMSE {runs[:5]}, median run {runs[5]}')
    rows_out = []
    for name in ('T2', 'T3', 'T3r', 'T3m'):
        b = base13[name]
        rows_out.append(dict(target=name, sigma='D+S (August)', run=13, **{k: v for k, v in b.items() if k not in ('by_system', 'contrib')}))
        for run in runs:
            e = res[run][name]
            rows_out.append(dict(target=name, sigma='D+S+K', run=run, **{k: v for k, v in e.items() if k not in ('by_system', 'contrib')}))
    O = pd.DataFrame(rows_out)
    pd.set_option('display.float_format', lambda v: f'{v:.3f}')
    print(O.to_string(index=False))
    print('\nT3r by system, design 13:'); print(res[13]['T3r']['by_system'].round(3).to_string())
    print('\nT3m by system, design 13:'); print(res[13]['T3m']['by_system'].round(3).to_string())
    print('\nT2 by system, design 13 (D+S+K = D+S: the term is centred):'); print(res[13]['T2']['by_system'].round(3).to_string())
    # sensitivity at design 13 (T3r)
    sens = []
    for lab, ch in [('final', {}), ('depth length x0.5', dict(l0=FINAL['l0'] * 0.5, l1=FINAL['l1'] * 0.5, lb0=FINAL['lb0'] * 0.5, lb1=FINAL['lb1'] * 0.5)),
                    ('depth length x2', dict(l0=FINAL['l0'] * 2, l1=FINAL['l1'] * 2, lb0=FINAL['lb0'] * 2, lb1=FINAL['lb1'] * 2)),
                    ('exponential depth corr', dict(nu=0.5, nub=0.5)), ('drift length 20 m', dict(ld=20.0)), ('drift length 200 m', dict(ld=200.0)),
                    ('no drift', dict(beta=0.0)), ('no between-borehole part', dict(u0=0.0)), ('no high-firn factor', dict(hf=1.0)),
                    ('not centred', dict(centre=False)), ('periodic share 0', dict(w=0.0)), ('periodic share 0.5', dict(w=0.5))]:
        p = dict(FINAL, **ch)
        e = evaluate(p, 13, names=('T3r', 'T2'))
        sens.append(dict(variant=lab, T3r=e['T3r']['ratio'], T3r_level=e['T3r']['level'], T3r_contrast=e['T3r']['contrast'],
                         T3r_nospike=e['T3r'].get('ratio_nospike'), T2=e['T2']['ratio']))
    print('\nsensitivity (design 13):'); print(pd.DataFrame(sens).to_string(index=False))
    pickle.dump(dict(O=O, res13=res[13], sens=pd.DataFrame(sens), runs=runs), open(W + 'final_ratios.pkl', 'wb'))
