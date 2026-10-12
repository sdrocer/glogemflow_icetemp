"""T1 and T2 posteriors on alps_s2 (as smoke_calibration.py), samples saved for the chapter-7 figure."""
import json
from pathlib import Path

import dill
import numpy as np

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']
priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
prior_sd = np.array([priors.refreeze_frac.std(), priors.insul_scale.std()])
save, summary = {}, {}
for name, target, offset in (('T2', 'depth', None), ('T1', 'entity', False)):
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=priors, track='depth', target=target,
                             offset_in_fit=offset, fixed_params={'advection_scale': 1.0}, fit_exclude=excl)
    sampler, flat = cal.run_mcmc(n_walkers=32, n_steps=3000, burn_in=1000, thin=5)
    save[name] = flat[:, :2]
    sd = flat[:, :2].std(axis=0)
    summary[name] = dict(mean=flat[:, :2].mean(axis=0).tolist(), sd=sd.tolist(), spread=(sd / prior_sd).tolist(),
                         corr=float(np.corrcoef(flat[:, 0], flat[:, 1])[0, 1]),
                         misfit_ratio=cal.mahalanobis_per_dof(flat.mean(axis=0)))
    print(name, summary[name], flush=True)
np.savez(OUT / 'posterior_T1_T2.npz', design=fx['theta'], **save)
(OUT / 'posterior_T1_T2.json').write_text(json.dumps(summary, indent=1))
