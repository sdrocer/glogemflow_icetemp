"""New code, target='entity' (default): identical to today's at 5 parameter sets.

Run: PYTHONPATH=<repo>/src python verify_entity.py
"""
import json
import sys
from pathlib import Path

import dill
import numpy as np

import icetemp.calibration.calibrator as calmod
from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.priors import Priors

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
print('code:', calmod.__file__)
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
ref = json.loads((OUT / 'entity_reference.json').read_text())
arr = np.load(OUT / 'entity_reference.npz')

cal = BayesianCalibrator(emulator=fx['emulator'], glaciers=fx['glaciers'], priors=Priors(), track='depth')
cal.fit_discrepancy()
ok = np.array_equal(cal._K_delta, arr['K_delta'])
for k, th in enumerate(arr['thetas']):
    r, S = cal.residual_and_covariance(th)
    same = (np.array_equal(r, arr[f'r{k}']) and np.array_equal(S, arr[f'S{k}'])
            and cal.log_likelihood(th) == ref['loglik'][k]
            and cal.mahalanobis_per_dof(th) == ref['mahalanobis'][k]
            and cal.log_posterior(th) == ref['logpost'][k])
    print(k, th, 'identical' if same else 'DIFFERENT', cal.log_likelihood(th), ref['loglik'][k])
    ok &= same
print('ALL IDENTICAL' if ok else 'MISMATCH')
sys.exit(0 if ok else 1)
