"""Step 5b fixture: emulator on the alps_s2 ensemble (29 runs with f_ref >= 0.1), built with
TODAY's package code (old_src), plus today's entity-target likelihood at 5 parameter sets.

Run: PYTHONPATH=old_src python build_fixture.py
"""
import json
import pickle
import time
from pathlib import Path

import dill
import numpy as np
import pandas as pd

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import Priors, PARAM_NAMES

REPO = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp')
TREES = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees')
OUT = TREES / 'step5b'
RESCORED = TREES / 'rescored' / 'alps_s2_step5'

THETAS = np.array([[0.5245, 1.0893, 1.054977],    # alps design 13
                   [0.30, 0.70, 1.00],
                   [0.80, 1.40, 1.00],
                   [0.55, 1.00, 0.80],
                   [0.20, 1.80, 1.30]])


def main():
    with open(TREES / 'step5a' / 'cache' / 'datahandler.pkl', 'rb') as f:
        data = pickle.load(f)
    known = {g.glacier_name: g for g in data['glaciers']}
    design = pd.read_csv(REPO / 'data' / 'bayescal' / 'alps_s2_w0' / 'design.csv')[list(PARAM_NAMES)]
    keep = np.flatnonzero(design.refreeze_frac.to_numpy() >= 0.1)
    outputs = []
    for i in keep:
        with open(RESCORED / f'design_{i:04d}.pkl', 'rb') as f:
            outputs.append(pickle.load(f))
    names = sorted(set().union(*outputs) & set(known), key=lambda n: (known[n].base_glacier_name,
                                                                         known[n].elevation))
    glaciers = [known[n] for n in names]
    theta = design.to_numpy()[keep]
    print(f'{len(glaciers)} entities, {len(keep)} runs')

    t0 = time.time()
    emu = Emulator(glaciers=glaciers, explained_variance=0.999)
    emu.fit(theta, outputs)
    emu.calibrate_variance(theta, outputs)
    print(f'emulator p={emu.p}, rows={len(emu.x_index_used())}, {time.time() - t0:.0f} s, '
          f'max_trusted_distance={emu.max_trusted_distance_:.3f}')
    with open(OUT / 'fixture.dill', 'wb') as f:
        dill.dump(dict(glaciers=glaciers, theta=theta, design_index=keep, outputs=outputs,
                       emulator=emu, systems=data['systems']), f)

    cal = BayesianCalibrator(emulator=emu, glaciers=glaciers, priors=Priors(), track='depth')
    cal.fit_discrepancy()
    ref = {'loglik': [], 'mahalanobis': [], 'logpost': []}
    arrays = {}
    for k, th in enumerate(THETAS):
        r, S = cal.residual_and_covariance(th)
        arrays[f'r{k}'], arrays[f'S{k}'] = r, S
        ref['loglik'].append(float(cal.log_likelihood(th)))
        ref['mahalanobis'].append(float(cal.mahalanobis_per_dof(th)))
        ref['logpost'].append(float(cal.log_posterior(th)))
    np.savez(OUT / 'entity_reference.npz', thetas=THETAS, K_delta=cal._K_delta, **arrays)
    (OUT / 'entity_reference.json').write_text(json.dumps(ref, indent=1))
    print(json.dumps(ref, indent=1))


if __name__ == '__main__':
    main()
