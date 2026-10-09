"""(c) Readings at the melting point. A point (T2) or row (T3r) is "at the melting point" when its
observed temperature is within 0.1 C of the pressure-melting point of its depth or above it:
obs >= PMP(z) - 0.1, PMP(z) = -0.000668 C/m x z (GloGEM's clamp, thermal_structure.py).
One-sided treatment: such a row is a lower bound. The model cannot be warmer than PMP; a model at
or above the reading (capped at PMP) costs nothing, a colder model is penalised as usual:
r = min(model - min(obs, PMP), 0), inside the same GLS quadratic form (rows keep their error
covariance). Compared: ordinary (baseline), one-sided, and the rows left out. Writes c05_melting.log/.pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
pd.set_option('display.width', 280, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
c2, c3 = ctx['T2'], ctx['T3r']
near = {}
for c in (c2, c3):
    p = c.cap
    d = c.y - p
    nr = d >= -0.1
    near[c.name] = nr
    print(f'{c.name}: rows {len(c.y)}; at the melting point (obs >= PMP - 0.1 C) {nr.sum()} ({nr.mean():.1%}); above PMP {int((d > 0).sum())} '
          f'(max +{d.max():.3f} C); depth {c.rows.depth[nr].min():.1f}-{c.rows.depth[nr].max():.1f} m')
    g = c.rows[nr].groupby('grp').agg(n=('depth', 'size'), zmin=('depth', 'min'), zmax=('depth', 'max'))
    g['obs_mean'] = c.rows[nr].groupby('grp').obs.mean() if 'obs' in c.rows else np.nan
    # model at those rows: design 13 and fraction of the 29 runs at PMP (within 0.01 C)
    Yk = c.Y[c.d.keep][:, nr]
    atp = (Yk >= c.cap[nr][None, :] - 0.01)
    g['model_at_pmp_d13'] = pd.Series(c.Y[13, nr] >= c.cap[nr] - 0.01, index=c.rows.index[nr]).groupby(c.rows.grp[nr]).sum()
    g['at_pmp_all29'] = pd.Series(atp.all(0), index=c.rows.index[nr]).groupby(c.rows.grp[nr]).sum()
    g['resid_d13'] = pd.Series(c.Y[13, nr] - c.y[nr], index=c.rows.index[nr]).groupby(c.rows.grp[nr]).mean()
    g['range_over_runs'] = pd.Series(np.ptp(Yk, axis=0), index=c.rows.index[nr]).groupby(c.rows.grp[nr]).median()
    print(g.round(3).to_string())
    allm = (c.Y[c.d.keep] >= c.cap[None, :] - 0.01)
    print(f'   all rows: model at PMP (within 0.01 C) at design 13: {int((c.Y[13] >= c.cap - 0.01).sum())}; in all 29 runs: {int(allm.all(0).sum())}; '
          f'of the near-PMP rows the model is colder than the reading by > 0.5 C at design 13: {int(((c.Y[13] - np.minimum(c.y, c.cap)) < -0.5)[nr].sum())}')
masks2 = {'T2 all': np.ones(len(c2.y), bool), 'T2 at melting point': near['T2'], 'T2 other': ~near['T2']}
variants = []
for c in (c2, c3):
    nr = near[c.name]
    variants += [('ordinary', c.name, np.ones(len(c.y), bool), None), ('one-sided', c.name, np.ones(len(c.y), bool), nr),
                 ('left out', c.name, ~nr, None)]
recs, fits = C.run_variants(ctx, variants, base, masks2)
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref',
                   'spread_f_ins', 'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd', 'd_f_ref_base_sd',
                   'd_f_ins_base_sd', 'mahal_base', 'chi2_n', 'rmse T2 all', 'rmse T2 at melting point', 'rmse T2 other'])
print(T.round(3).to_string(index=False))
# at each variant's optimum: how many near-PMP rows are penalised (model colder than min(obs, PMP))
for (lab, tname), r in fits.items():
    c = ctx[tname]
    nr = near[tname]
    for kind, s in r.items():
        th = np.array([s['map_f_ref'], s['map_f_ins'], 1.0])
        G = c.sur.predict(th)[0] if kind == 'surface' else c.emu.predict(th)[0][0]
        yb = np.minimum(c.y, c.cap)
        pen = (G < yb - 1e-9) & nr
        print(f'{tname} {lab:9s} {kind:8s}: near-PMP rows with the model colder than the reading: {pen.sum()} of {nr.sum()}, '
              f'their mean shortfall {np.mean((G - yb)[pen]):+.3f} C; model at PMP (within 0.01): {int(((G >= c.cap - 0.01) & nr).sum())}')
pickle.dump(dict(table=T, fits=fits, near={k: v for k, v in near.items()}), open(C.OUT / 'c05_melting.pkl', 'wb'))
