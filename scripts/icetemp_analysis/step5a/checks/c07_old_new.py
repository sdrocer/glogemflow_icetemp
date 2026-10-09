"""(g) Old against new readings: points/rows measured before 1979 (all their readings before
1 Jan 1979) against later ones. Residuals (bias, RMSE, misfit ratio r^T Sigma^-1 r / n under the
target's error model) at design 13 and at the T2 / T3r optimum; the pairs where one site was
measured in both periods (Col du Dome: Bossons 4275-4285 m in 1973-76 vs 4245-4305 m in
1994-2011; Colle Gnifetti: Grenz@4455 in 1977 vs later); the fit without the old rows; and
option C (28 Sep): an extra model-error variance on the old rows, sized so that their misfit
ratio equals the new rows'. Writes c07_old_new.log / .pkl."""
import sys, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
from scipy.optimize import brentq
import c3lib as C
ta = C.ta
pd.set_option('display.width', 280, 'display.max_columns', 40, 'display.max_rows', 200)
ctx = C.context()
base = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
c2, c3 = ctx['T2'], ctx['T3r']
CUT = 197901
old2 = (c2.rows.ym_last < CUT).to_numpy(); new2 = (c2.rows.ym_first >= CUT).to_numpy()
old3 = (c3.rows.ym_last < CUT).to_numpy(); new3 = (c3.rows.ym_first >= CUT).to_numpy()
print(f'T2: old points {old2.sum()}, new {new2.sum()}, mixed {(~old2 & ~new2).sum()}; T3r: old rows {old3.sum()}, new {new3.sum()}, mixed {(~old3 & ~new3).sum()}')
print('T2 old points by entity:', c2.rows[old2].entity.value_counts().sort_index().to_dict())
print('T2 mixed points:', c2.rows[~old2 & ~new2][['entity', 'depth', 'ym_first', 'ym_last']].to_string())
print('T3r old rows by entity:', c3.rows[old3].entity.value_counts().sort_index().to_dict())
HEF = 'Hintereisferner tongue, ice (< 3000 m)'


def misfit(c, r, m):
    S = c.Sigma[np.ix_(m, m)]
    L = np.linalg.cholesky(S)
    w = np.linalg.solve(L, r[m])
    return float(w @ w / m.sum())


def stats(c, G, m):
    r = G - c.y
    return dict(n=int(m.sum()), bias=float(np.mean(r[m])), rmse=C.rmse(G[m], c.y[m]), ratio=misfit(c, r, m))


rows = []
for c, old, new in ((c2, old2, new2), (c3, old3, new3)):
    hef = (c.rows.grp == HEF).to_numpy()
    cdd_old = c.rows.entity.isin(['Glacier des Bossons@4275', 'Glacier des Bossons@4285']).to_numpy() & old
    cdd_new = c.rows.entity.isin(['Glacier des Bossons@4245', 'Glacier des Bossons@4255', 'Glacier des Bossons@4295',
                                  'Glacier des Bossons@4305', 'Glacier des Bossons@4265']).to_numpy() & new
    cg_old = (c.rows.entity == 'Grenzgletscher@4455').to_numpy() & old
    cg_new = (c.rows.entity == 'Grenzgletscher@4455').to_numpy() & new
    sets = [('old (before 1979)', old), ('new (1979 and later)', new), ('old without Hintereisferner tongue', old & ~hef),
            ('Hintereisferner tongue (all old)', hef), ('Col du Dome old (Bossons 4275-4285, 1973-76)', cdd_old),
            ('Col du Dome new (Bossons 4245-4305, 1994-2011)', cdd_new), ('Colle Gnifetti 4455 old (1977)', cg_old),
            ('Colle Gnifetti 4455 new (1990-2021)', cg_new)]
    for kind in ('design 13', 'surface optimum', 'emulator optimum'):
        if kind == 'design 13':
            G = c.Y[13]
        else:
            s = base[c.name]['surface' if kind.startswith('surface') else 'emulator']
            th = np.array([s['map_f_ref'], s['map_f_ins'], 1.0])
            G = c.sur.predict(th)[0] if kind.startswith('surface') else c.emu.predict(th)[0][0]
        for lab, m in sets:
            if m.sum() == 0:
                continue
            rows.append(dict(target=c.name, at=kind, set=lab, **stats(c, G, m)))
R = pd.DataFrame(rows)
print(R.round(3).to_string(index=False))
# per old entity at the T2 surface optimum
s = base['T2']['surface']; G = c2.sur.predict(np.array([s['map_f_ref'], s['map_f_ins'], 1.0]))[0]
print('\nT2 per entity with old points (surface optimum):')
print(pd.DataFrame([dict(entity=e, **stats(c2, G, (c2.rows.entity == e).to_numpy() & old2)) for e in sorted(c2.rows[old2].entity.unique())]).round(3).to_string(index=False))
for e in ('Glacier des Bossons@4275', 'Glacier des Bossons@4285'):
    m = (c2.rows.entity == e).to_numpy()
    print(e, 'points (depth, obs, model at optimum, first-last month):', [(round(z, 1), round(o, 2), round(g, 2), a, b) for z, o, g, a, b in zip(c2.rows.depth[m], c2.y[m], G[m], c2.rows.ym_first[m], c2.rows.ym_last[m])])
# fit without the old rows, and option C
base_r = pickle.load(open(C.OUT / 'baseline.pkl', 'rb'))
variants = [('without old rows', 'T2', ~old2, None), ('without old rows', 'T3r', ~old3, None),
            ('without old rows except Hintereisferner tongue', 'T2', ~(old2 & ~(c2.rows.grp == HEF).to_numpy()), None)]
recs, fits = C.run_variants(ctx, variants, base_r, {'T2 all': np.ones(len(c2.y), bool), 'T2 old': old2, 'T2 new': new2})
# option C: extra variance v on the old rows, sized at the surface optimum so that the old rows' misfit ratio
# equals the new rows'; as an offset per old entity (correlated within it) and, if an offset cannot reach it,
# as an independent variance per row
for c, old, new in ((c2, old2, new2), (c3, old3, new3)):
    s = base[c.name]['surface']
    G = c.sur.predict(np.array([s['map_f_ref'], s['map_f_ins'], 1.0]))[0]
    r = G - c.y
    target = misfit(c, r, new)
    e = c.rows.e.to_numpy()
    same = ((e[:, None] == e[None, :]) & old[:, None] & old[None, :]).astype(float)
    diag = np.diag(old.astype(float))

    def ratio_old(v, K):
        S = (c.Sigma + v * K)[np.ix_(old, old)]
        w = np.linalg.solve(np.linalg.cholesky(S), r[old])
        return float(w @ w / old.sum())
    r0 = ratio_old(0.0, same)
    print(f'\n{c.name}: misfit ratio old {r0:.2f} vs new {target:.2f} at the surface optimum; with an offset variance of 20 / 100 C^2 on old '
          f'entities: {ratio_old(20.0, same):.2f} / {ratio_old(100.0, same):.2f}; with an independent row variance of 2 / 10 C^2: '
          f'{ratio_old(2.0, diag):.2f} / {ratio_old(10.0, diag):.2f}')
    for lab, K in (('offset per old entity', same), ('independent per old row', diag)):
        if r0 <= target:
            print(f'   {lab}: not needed (old rows not worse)'); continue
        if ratio_old(500.0, K) > target:
            print(f'   {lab}: cannot reach the new rows\' ratio (floor {ratio_old(500.0, K):.2f} at 500 C^2)'); continue
        v = brentq(lambda x: ratio_old(x, K) - target, 0.0, 500.0)
        Sg = c.Sigma + v * K
        r_ = C.fit_rows(c, np.ones(len(c.y), bool), Sigma=Sg)
        fits[(f'option C, {lab}', c.name)] = C.strip(r_)
        rr = C.row(f'option C: +{v:.2f} C^2 {lab}', c.name, r_, base_r[c.name])
        for rec in rr:
            rm = C.rmse_block(c2, r_[rec['method']], {'T2 all': np.ones(len(c2.y), bool), 'T2 old': old2, 'T2 new': new2}, rec['method'])
            rec.update({f'rmse {k}': vv for k, vv in rm.items()})
        recs += rr
        print(f'   {lab}: {v:.2f} C^2')
T = C.table(recs, ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr', 'spread_f_ref',
                   'spread_f_ins', 'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd', 'd_f_ref_base_sd',
                   'd_f_ins_base_sd', 'mahal_base', 'chi2_n', 'rmse T2 all', 'rmse T2 old', 'rmse T2 new'])
print(T.round(3).to_string(index=False))
pickle.dump(dict(residuals=R, table=T, fits=fits), open(C.OUT / 'c07_old_new.pkl', 'wb'))
