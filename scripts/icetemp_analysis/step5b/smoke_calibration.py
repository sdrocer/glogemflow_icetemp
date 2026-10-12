"""Step 5b smoke calibration on alps_s2 (7 glaciers, 66 entities, 29 runs, 3-parameter emulator at
99.9 %): T2 target, f_adv fixed at 1.0, f_ref prior U(0.1, 1), Colle Gnifetti below 40 m out of the
fit. Grid, mode, MCMC, gate, held-out test by glacier system; T1 (entity, no offset) for comparison.
"""
import json
import time
from pathlib import Path

import dill
import numpy as np

from icetemp.calibration import adoption
from icetemp.calibration.baselines import grid_search_all
from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors
from icetemp.calibration.validation import Validator

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']
priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
FIXED = {'advection_scale': 1.0}
prior_sd = np.array([priors.refreeze_frac.std(), priors.insul_scale.std()])
res = {}


def make(target, offset=None):
    return BayesianCalibrator(emulator=emu, glaciers=gl, priors=priors, track='depth', target=target,
                              offset_in_fit=offset, fixed_params=FIXED, fit_exclude=excl)


for name, target, offset in (('T2', 'depth', None), ('T1', 'entity', False)):
    cal = make(target, offset)
    t0 = time.time()
    g_theta, g_lp = cal.grid_search()
    m_theta = cal.find_map()
    sampler, flat = cal.run_mcmc(n_walkers=32, n_steps=3000, burn_in=1000, thin=5)
    summ, tau = cal.convergence_diagnostics(sampler, 1000, 5)
    mean = flat.mean(axis=0)
    mcmc_lp = max(cal.log_posterior(t) for t in flat[::max(len(flat) // 2000, 1)])
    sd = flat[:, :2].std(axis=0)
    corr = float(np.corrcoef(flat[:, 0], flat[:, 1])[0, 1])
    maha = cal.mahalanobis_per_dof(mean)
    gate = adoption.check_gate(emu, mean, r_hat_max=float(summ['r_hat'].max()), mahalanobis=maha,
                               mode_gap=float(g_lp - mcmc_lp))
    res[name] = dict(rows=len(cal._row_obs), grid=g_theta.tolist(), grid_lp=g_lp, map=m_theta.tolist(),
                     post_mean=mean.tolist(), post_sd=sd.tolist(), spread=(sd / prior_sd).tolist(),
                     corr=corr, r_hat=float(summ['r_hat'].max()), ess=float(summ['ess_bulk'].min()),
                     mode_gap=float(g_lp - mcmc_lp), misfit_ratio=maha, n_nearby=gate.n_nearby,
                     gate=gate.passed, gate_reasons=gate.reasons, seconds=round(time.time() - t0))
    print(name, json.dumps(res[name], indent=1), flush=True)

cal = make('depth')
val = Validator(glaciers=gl, calib_df=grid_search_all(gl), calibrator=cal)
t0 = time.time()
points, folds = val.run_parameters()
text, pooled, ok = val.parameters_summary(points, folds, margin=0.10)
print(text, f'({time.time() - t0:.0f} s)')
print(folds[['system', 'n_points', 'ko_refreeze_frac', 'ko_insul_scale', 't1_refreeze_frac',
             't1_insul_scale', 'ko_in_design', 't1_in_design']].round(3).to_string(index=False))
points.to_csv(OUT / 'smoke_heldout_points.csv', index=False)
folds.to_csv(OUT / 'smoke_heldout_folds.csv', index=False)
res['heldout'] = dict(pooled=pooled, rule=ok)
(OUT / 'smoke_calibration.json').write_text(json.dumps(res, indent=1, default=float))

# the old comparison (baselines with the offset, classification) still runs with the depth target
t0 = time.time()
results = val.run()
print(val.summary(results)[0], f'({time.time() - t0:.0f} s)')
