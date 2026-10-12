"""target='depth' on the alps_s2 fixture: fast (Woodbury) vs full likelihood, speed, misfit ratio
at design 13 (5a: 0.605 from the real run, Sigma = D + S), with and without Colle Gnifetti > 40 m.
"""
import time
from pathlib import Path

import dill
import numpy as np

from icetemp.calibration.calibrator import BayesianCalibrator, _gaussian_logpdf
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
emu, gl = fx['emulator'], fx['glaciers']
j13 = int(np.flatnonzero(fx['design_index'] == 13)[0])
th13 = fx['theta'][j13]
excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
print(f'design 13 = {th13}; Colle Gnifetti > 40 m: {len(excl)} points')

for label, fe in (('all points', None), ('without CG > 40 m', excl)):
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth',
                             fixed_params={'advection_scale': 1.0}, fit_exclude=fe)
    n = len(cal._row_obs)
    # real run 13, Sigma = D + S (5a's misfit ratio)
    x_used = emu.x_index_used()
    real = []
    for row in cal._row_emu_idx:
        g, d = x_used[row]
        depths, T = fx['outputs'][j13][g]
        real.append(T[int(np.argmin(np.abs(np.asarray(depths) - d)))])
    r = cal._row_obs - np.array(real)
    DS = np.diag(cal._row_sigma_obs2) + cal._row_model_covariance()
    ratio_real = r @ np.linalg.solve(DS, r) / n
    print(f'\n{label}: {n} rows, {len(cal._calib_glacier_names)} entities')
    print(f'  misfit ratio at run 13, real output, D+S: {ratio_real:.3f}')
    print(f'  misfit ratio at design 13, emulator, D+S+trunc+E: {cal.mahalanobis_per_dof(th13):.3f} '
          f'(truncation variance {emu.truncation_variance():.4f})')
    for th in (th13, [0.3, 0.7, 1.0], [0.8, 1.4, 1.0], [0.55, 1.0, 1.0], [0.2, 1.8, 1.3]):
        fast = cal.log_likelihood(th)
        full = _gaussian_logpdf(*cal.residual_and_covariance(th))
        print(f'  theta {np.round(th, 3)}: fast {fast:.10f} full {full:.10f} '
              f'diff {fast - full:.2e} rel {abs(fast - full) / abs(full):.1e}')
    t0 = time.perf_counter()
    for _ in range(300):
        cal.log_likelihood(th13)
    t_fast = (time.perf_counter() - t0) / 300
    t0 = time.perf_counter()
    for _ in range(20):
        _gaussian_logpdf(*cal.residual_and_covariance(th13))
    t_full = (time.perf_counter() - t0) / 20
    t0 = time.perf_counter()
    for _ in range(300):
        emu.predict_factor(th13)
    t_emu = (time.perf_counter() - t0) / 300
    print(f'  per evaluation: fast {1e3 * t_fast:.2f} ms (emulator alone {1e3 * t_emu:.2f} ms), '
          f'full {1e3 * t_full:.1f} ms')
