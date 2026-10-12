"""Log-posterior grids on alps_s2 (f_adv 1.0) for the chapter-7 figure: T1 (entity, no offset), T2 on
the emulator mean only, T2 with the emulator uncertainty varying; trust region and design runs."""
from pathlib import Path
import dill, numpy as np
from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']
priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])


class NoE:
    def __init__(self, e):
        self.e = e
    def __getattr__(self, k):
        return getattr(self.e, k)
    def predict_factor(self, th):
        m, L = self.e.predict_factor(th)
        return m, np.zeros_like(L)


fr, fi = np.linspace(0.1, 1.0, 91), np.linspace(0.3, 2.0, 86)
out = dict(fr=fr, fi=fi, design=fx['theta'], trusted=emu.max_trusted_distance_,
           lo=emu.theta_train_.min(axis=0), hi=emu.theta_train_.max(axis=0))
out['nn'] = np.array([[emu._nn_distance(np.array([[a, b, 1.0]]))[0] for b in fi] for a in fr])
for name, target, offset, e in (('T1', 'entity', False, emu), ('T2_mean', 'depth', None, NoE(emu)),
                                ('T2_E', 'depth', None, emu)):
    cal = BayesianCalibrator(emulator=e, glaciers=gl, priors=priors, track='depth', target=target,
                             offset_in_fit=offset, fixed_params={'advection_scale': 1.0}, fit_exclude=excl)
    out[name] = np.array([[priors.logpdf((a, b, 1.0)) + cal.log_likelihood(np.array([a, b, 1.0]))
                           for b in fi] for a in fr])
    print(name, 'done', flush=True)
np.savez(OUT / 'surface_grids.npz', **out)
