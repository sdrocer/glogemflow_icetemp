"""Held-out test on alps_s2 with the fits on the emulator mean only (E = 0), as smoke_calibration.py."""
from pathlib import Path
import dill, numpy as np
from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors
from icetemp.calibration.validation import Validator

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']


class NoE:
    def __init__(self, e):
        self.e = e
    def __getattr__(self, k):
        return getattr(self.e, k)
    def predict_factor(self, th):
        m, L = self.e.predict_factor(th)
        return m, np.zeros_like(L)


cal = BayesianCalibrator(emulator=NoE(emu), glaciers=gl, priors=Priors(bounds_override={'refreeze_frac': [0.1, 1.0]}),
                         track='depth', target='depth', fixed_params={'advection_scale': 1.0},
                         fit_exclude=fit_exclusions(gl, ['colle_gnifetti_below_40m']))
val = Validator(glaciers=gl, calib_df=None, calibrator=cal)
points, folds = val.run_parameters()
text, pooled, ok = val.parameters_summary(points, folds, margin=0.10)
print(text)
print(folds[['system', 'n_points', 'ko_refreeze_frac', 'ko_insul_scale', 't1_refreeze_frac', 't1_insul_scale',
             'ko_nn_dist']].round(3).to_string(index=False))
points.to_csv(OUT / 'heldout_noE_points.csv', index=False)
folds.to_csv(OUT / 'heldout_noE_folds.csv', index=False)
