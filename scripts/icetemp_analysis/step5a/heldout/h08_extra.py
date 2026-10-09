"""Extra checks: (1) how well the surrogates predict the T2 RMSE of real runs near the T2/T3 optima
(leave-one-run-out); (2) per-system T2 RMSE at each target's in-sample fit; (3) mean misfit
(model - obs) by system and depth band at the level-only and shape-only optima; (4) consistency
distance with the level's covariance inflated by its own misfit ratio (Birge)."""
import sys, pickle, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

pd.set_option('display.width', 250, 'display.max_columns', 40)
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
t2 = T['T2']
keep = design['keep'].to_numpy()
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
idx = np.flatnonzero(keep)
loo = pickle.load(open(H.OUT / 'h01_loo.pkl', 'rb'))['loos']['T2']
true_rmse = np.sqrt(((t2.Y[idx] - t2.y) ** 2).mean(1))
print('== (1) leave-one-run-out: surrogate-predicted vs true T2 RMSE of the real runs [C] ==')
rows = []
for opt_name, opt in (('T1 fit', (0.53, 1.005)), ('T2 fit', (0.34, 0.73)), ('T3 fit', (0.26, 0.645))):
    dist = np.sqrt((((th3[idx, :2] - opt) / H.SD) ** 2).sum(1))
    for j in np.argsort(dist)[:3]:
        r = dict(near=opt_name, run=int(idx[j]), f_ref=th3[idx[j], 0], f_ins=th3[idx[j], 1], f_adv=th3[idx[j], 2],
                 dist_sd=dist[j], true=true_rmse[j])
        for k in ('C2Q', 'Q3', 'EMU999'):
            r[k] = float(np.sqrt(((t2.Y[idx[j]] + loo[k][j] - t2.y) ** 2).mean()))
        rows.append(r)
X = pd.DataFrame(rows)
print(X.round(3).to_string(index=False))
allerr = {k: np.sqrt(((t2.Y[idx] + loo[k] - t2.y) ** 2).mean(1)) - true_rmse for k in ('C2Q', 'Q3', 'EMU999')}
print('all 29 runs, predicted minus true T2 RMSE: ' + ', '.join(
    f'{k} median {np.median(v):+.3f}, mean |d| {np.abs(v).mean():.3f}, max |d| {np.abs(v).max():.3f}' for k, v in allerr.items()))

# (2) per-system RMSE at the in-sample fits (response surface)
rs = H.RespSurf('C2Q').fit(th3[keep], t2.Y[keep])
fits = {'T1 fit (0.530, 1.005)': (0.53, 1.005), 'T2 fit (0.340, 0.730)': (0.34, 0.73), 'T3 fit (0.260, 0.645)': (0.26, 0.645),
        'level only (0.530, 1.005)': (0.53, 1.005), 'T2 shape only (0.240, 0.845)': (0.24, 0.845),
        'T3 shape only (0.200, 0.830)': (0.20, 0.83)}
pred = {k: rs.predict(np.array([v]))[0] for k, v in fits.items()}
pred['run 13 (real)'] = t2.Y[13]
zb = H.ta.depth_band(t2.rows['depth'].to_numpy(float))
sysn = t2.system
print('\n== (2) T2 RMSE by system at the in-sample fits (response surface) [C] ==')
tab = {}
for k in ('run 13 (real)', 'T1 fit (0.530, 1.005)', 'T2 fit (0.340, 0.730)', 'T3 fit (0.260, 0.645)'):
    r = pred[k] - t2.y
    tab[k] = {s: float(np.sqrt((r[sysn == s] ** 2).mean())) for s in sorted(set(sysn))}
    tab[k]['all'] = float(np.sqrt((r ** 2).mean()))
tab = pd.DataFrame(tab)
tab['points'] = pd.Series({s: int((sysn == s).sum()) for s in sorted(set(sysn))} | {'all': t2.n})
print(tab.round(3).to_string())
print('\n== (3) mean misfit model - obs [C] by system and depth band, at the level-only and shape-only optima ==')
lab = np.where((sysn == 'Monte Rosa') & (t2.rows['band'].to_numpy() >= 3900), 'Monte Rosa firn >= 3900 m',
               np.where(sysn == 'Monte Rosa', 'Monte Rosa < 3900 m', sysn))
rec = []
for k in ('level only (0.530, 1.005)', 'T2 shape only (0.240, 0.845)', 'T3 shape only (0.200, 0.830)'):
    r = pred[k] - t2.y
    for g in sorted(set(lab)):
        for b in ('0-5 m', '5-15 m', '15-30 m', '> 30 m'):
            m = (lab == g) & (zb == b)
            if m.sum():
                rec.append(dict(at=k, group=g, band=b, n=int(m.sum()), bias=float(r[m].mean())))
B = pd.DataFrame(rec)
P = B.pivot_table(index=['group', 'band'], columns='at', values='bias').round(2)
P['n'] = B.groupby(['group', 'band']).n.first()
print(P.to_string())
# (4) Birge-scaled consistency
C6 = json.load(open(H.OUT / 'h06_consistency.json'))
H7 = json.load(open(H.OUT / 'h07_honest.json'))
R7 = pd.read_pickle(H.OUT / 'h07_honest.pkl')
t1_ratio = float(R7[(R7.target == 'T1') & (R7.run == H7['best']['T1'])].ratio.iloc[0])


def cov_sd(d):
    return np.array([[d['sd_ref'] ** 2, d['corr'] * d['sd_ref'] * d['sd_ins']],
                     [d['corr'] * d['sd_ref'] * d['sd_ins'], d['sd_ins'] ** 2]])


print(f'\n== (4) consistency with the level covariance inflated by the T1 misfit ratio at its best run ({t1_ratio:.3f}) ==')
out4 = {}
for key in ('RS/E', 'EMU/E'):
    mm = C6['main'][key]
    L = mm['level (T1)']
    for tg in ('T2', 'T3'):
        Sh = mm[f'{tg} shape']
        d = (np.array(L['th']) - np.array(Sh['th'])) / H.SD
        for f in (1.0, t1_ratio, 4.0):
            Cs = f * cov_sd(L) + cov_sd(Sh)
            out4[f'{key} {tg} x{f:.2f}'] = float(np.sqrt(d @ np.linalg.solve(Cs, d)))
print(json.dumps({k: round(v, 2) for k, v in out4.items()}, indent=1))
json.dump(dict(loo_near=X.to_dict('records'), rmse_by_system=tab.to_dict(), bias=B.to_dict('records'), birge=out4),
          open(H.OUT / 'h08_extra.json', 'w'), indent=1, default=float)
