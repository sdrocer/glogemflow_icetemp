"""How good is the RMSE below the seasonal layer? (question of 2026-10-12)

alps_s2 fixture (7 glaciers, 29 runs with f_ref >= 0.1, emulator at 99.9 %), current scheme: T2 on
the emulator mean, f_adv 1.0, Colle Gnifetti below 40 m out of the fit, held out one glacier system
at a time. "Deep" = point depth > 15 m (thermal_structure.SEASONAL_DEPTH_M); 10 and 20 m as checks.

A  the 29 real runs (no fit, no emulator): pooled RMSE and bias, all / shallow / deep points
B  held-out RMSE on the deep points: T2 and T1 fitted on all depths, T2 and T1 fitted on the deep
   points only, and the prior mean; per glacier system
C  in-sample fits on all systems: parameters and whether f_ref and f_ins are told apart
"""
import json
from pathlib import Path

import dill
import numpy as np
import pandas as pd

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors
from icetemp.calibration.thermal_structure import SEASONAL_DEPTH_M
from icetemp.calibration.validation import Validator

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']
priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
cg = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
FIXED = {'advection_scale': 1.0}
prior_sd = {'refreeze_frac': priors.refreeze_frac.std(), 'insul_scale': priors.insul_scale.std()}
res = {}


def rmse(x):
    return float(np.sqrt(np.mean(np.square(x)))) if len(x) else float('nan')


# ── A: real runs ──────────────────────────────────────────────────────────────────────────────
rows = [(g.glacier_name, float(d), float(t), (g.glacier_name, float(d)) in cg)
        for g in gl for d, t in zip(g.depths, g.T_obs)]
P = pd.DataFrame(rows, columns=['entity', 'depth', 'obs', 'cg_deep'])
x_used = set(emu.x_index_used())
P = P[[(e, d) in x_used for e, d in zip(P.entity, P.depth)]].reset_index(drop=True)
model = np.array([[out[e][1][int(np.argmin(np.abs(np.asarray(out[e][0]) - d)))] for e, d in zip(P.entity, P.depth)]
                  for out in fx['outputs']])                   # (29 runs, points)
R = model - P.obs.to_numpy()
j13 = int(np.flatnonzero(fx['design_index'] == 13)[0])
recA = []
for label, m in (('all points (fit set)', ~P.cg_deep),
                 (f'shallow, <= {SEASONAL_DEPTH_M:.0f} m', ~P.cg_deep & (P.depth <= SEASONAL_DEPTH_M)),
                 (f'deep, > {SEASONAL_DEPTH_M:.0f} m', ~P.cg_deep & (P.depth > SEASONAL_DEPTH_M)),
                 ('deep, > 10 m', ~P.cg_deep & (P.depth > 10)), ('deep, > 20 m', ~P.cg_deep & (P.depth > 20)),
                 (f'deep incl. Colle Gnifetti > 40 m', P.depth > SEASONAL_DEPTH_M),
                 ('Colle Gnifetti > 40 m only', P.cg_deep)):
    m = m.to_numpy()
    r = np.sqrt(np.mean(R[:, m] ** 2, axis=1))
    b = R[:, m].mean(axis=1)
    k = int(np.argmin(r))
    recA.append(dict(points=label, n=int(m.sum()), best_run=int(fx['design_index'][k]), rmse_best=r[k],
                     bias_best=b[k], rmse_run13=r[j13], bias_run13=b[j13], rmse_median=float(np.median(r))))
A = pd.DataFrame(recA)
print('== A: the 29 real runs (model minus boreholes) ==')
print(A.round(3).to_string(index=False))
res['A'] = recA

# ── B and C: calibrations ─────────────────────────────────────────────────────────────────────
deep_only = cg | {(g.glacier_name, float(d)) for g in gl for d in g.depths if d <= SEASONAL_DEPTH_M}


def make(target, fit_exclude, offset=None):
    return BayesianCalibrator(emulator=emu, glaciers=gl, priors=priors, track='depth', target=target,
                              offset_in_fit=offset, fixed_params=FIXED, fit_exclude=fit_exclude,
                              emulator_error='mean')


held = {}
for fit_name, excl in (('fitted on all depths', cg), ('fitted below 15 m only', deep_only)):
    val = Validator(glaciers=gl, calib_df=None, calibrator=make('depth', excl))
    points, folds = val.run_parameters()
    points['deep'] = points.depth > SEASONAL_DEPTH_M
    points['cg_deep'] = [(e, d) in cg for e, d in zip(points.entity, points.depth)]
    held[fit_name] = (points, folds)
    points.to_csv(OUT / f"deep_layer_heldout_{'all' if excl is cg else 'deep'}.csv", index=False)

print(f'\n== B: held-out RMSE / bias on the points below {SEASONAL_DEPTH_M:.0f} m (Colle Gnifetti > 40 m '
      f'apart), through the emulator, one glacier system per fold ==')
recB = []
for fit_name, (points, folds) in held.items():
    for scored, m in (('deep, fit set', points.deep & ~points.cg_deep), ('all fit points', ~points.cg_deep),
                      ('Colle Gnifetti > 40 m', points.cg_deep)):
        sub = points[m]
        for meth, label in (('ko', 'T2'), ('t1', 'T1'), ('prior', 'prior mean')):
            r = sub[f'pred_{meth}'] - sub.obs
            recB.append(dict(fit=fit_name, method=label, scored=scored, n=len(sub), rmse=rmse(r), bias=float(r.mean())))
B = pd.DataFrame(recB)
print(B.pivot_table(index=['scored', 'method'], columns='fit', values='rmse').round(3).to_string())
print('\nbias (model minus boreholes) [C]:')
print(B.pivot_table(index=['scored', 'method'], columns='fit', values='bias').round(3).to_string())
res['B'] = recB

print('\nper glacier system, deep points of the fit set (RMSE, degC):')
per = []
for fit_name, (points, _) in held.items():
    sub = points[points.deep & ~points.cg_deep]
    for system, df in sub.groupby('system'):
        rec = dict(fit=fit_name, system=system[:40], n=len(df))
        rec.update({lab: rmse(df[f'pred_{m}'] - df.obs) for m, lab in (('ko', 'T2'), ('t1', 'T1'), ('prior', 'prior'))})
        per.append(rec)
print(pd.DataFrame(per).round(2).to_string(index=False))
res['per_system'] = per

print('\nfold parameters (f_ref, f_ins):')
for fit_name, (_, folds) in held.items():
    print(f'  {fit_name}:')
    print(folds[['system', 'ko_refreeze_frac', 'ko_insul_scale', 't1_refreeze_frac', 't1_insul_scale']]
          .round(3).to_string(index=False))

print('\n== C: in-sample fits on all 6 systems (grid posterior; told apart: corr <= 0.8, f_ins spread <= 0.5) ==')
recC = []
for fit_name, excl in (('all depths', cg), ('below 15 m only', deep_only)):
    for target, offset, label in (('depth', None, 'T2'), ('entity', False, 'T1')):
        cal = make(target, excl, offset)
        g = cal.grid_moments()
        th = cal._expand_theta(g['mean'])
        mean = emu.predict(np.atleast_2d(th))[0][0]
        rows_idx = cal._row_emu_idx
        r_fit = mean[rows_idx] - cal._row_obs
        P_fit = [emu.x_index_used()[i][1] for i in rows_idx]
        recC.append(dict(fit=fit_name, target=label, n_rows=len(rows_idx), f_ref=g['mean'][0], f_ins=g['mean'][1],
                         spread_f_ref=g['sd'][0] / prior_sd['refreeze_frac'],
                         spread_f_ins=g['sd'][1] / prior_sd['insul_scale'], corr=g['corr'][0, 1],
                         told_apart=bool(abs(g['corr'][0, 1]) <= 0.8 and g['sd'][1] / prior_sd['insul_scale'] <= 0.5),
                         rmse_fit_rows=rmse(r_fit), misfit_ratio=cal.mahalanobis_per_dof(th)))
C = pd.DataFrame(recC)
print(C.round(3).to_string(index=False))
res['C'] = recC
(OUT / 'deep_layer.json').write_text(json.dumps(res, indent=1, default=float))

# ── D: in-sample RMSE below 15 m through the emulator, at each fit and at the prior mean ─────────
print(f'\n== D: in-sample RMSE / bias on the deep points of the fit set, through the emulator ==')
x_list = emu.x_index_used()
deep_rows = [i for i, (e, d) in enumerate(x_list)
             if d != 'basal' and d > SEASONAL_DEPTH_M and (e, d) not in cg]
obs = {(g.glacier_name, float(d)): float(t) for g in gl for d, t in zip(g.depths, g.T_obs)}
y = np.array([obs[x_list[i]] for i in deep_rows])
pm = priors.as_dict()
thetas = {'prior mean': np.array([pm['refreeze_frac'].mean(), pm['insul_scale'].mean(), 1.0]),
          'run 13 parameters, f_adv 1.0': np.array([0.5245, 1.0893, 1.0])}
for rec in recC:
    thetas[f"{rec['target']} fitted on {rec['fit']}"] = np.array([rec['f_ref'], rec['f_ins'], 1.0])
recD = []
for name, th in thetas.items():
    r = emu.predict(np.atleast_2d(th))[0][0][deep_rows] - y
    recD.append(dict(theta=name, f_ref=th[0], f_ins=th[1], n=len(y), rmse=rmse(r), bias=float(r.mean())))
r13 = emu.predict(np.atleast_2d(fx['theta'][j13]))[0][0][deep_rows] - y
recD.append(dict(theta='run 13 as run (f_adv 1.055), emulator', f_ref=fx['theta'][j13][0], f_ins=fx['theta'][j13][1],
                 n=len(y), rmse=rmse(r13), bias=float(r13.mean())))
print(pd.DataFrame(recD).round(3).to_string(index=False))
res['D'] = recD
(OUT / 'deep_layer.json').write_text(json.dumps(res, indent=1, default=float))
