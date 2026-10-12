"""Fixed parameters: the design spans the free ones only, design bounds follow the emulator's
columns, and a fixed parameter never moves in the search."""
import numpy as np
import pytest
from scipy.stats import qmc

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import PARAM_NAMES, Priors


def test_fixed_advection_never_moves(setup):
    gl, _, _, emu = setup
    fixed = {'advection_scale': 1.0}
    design = DesignSampler(Priors(bounds_override={'refreeze_frac': (0.1, 1.0)}), seed=1,
                           fixed_params=fixed).sample(16)
    assert np.all(design[:, 2] == 1.0) and design[:, 0].min() >= 0.1
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', fixed_params=fixed)
    cal.fit_discrepancy()
    assert cal._free_param_names() == ('refreeze_frac', 'insul_scale')
    assert cal.find_map(n_starts=4)[2] == 1.0
    assert cal.grid_search(n_per_dim=9)[0][2] == 1.0
    _, flat = cal.run_mcmc(n_walkers=8, n_steps=40, burn_in=10, thin=1)
    assert np.all(flat[:, 2] == 1.0) and flat[:, :2].std(axis=0).min() > 0


def test_design_without_fixed_params_unchanged():
    pri = Priors()
    u = qmc.LatinHypercube(d=3, scramble=True, seed=42).random(n=10)
    expect = np.column_stack([pri.as_dict()[n].ppf(u[:, k]) for k, n in enumerate(PARAM_NAMES)])
    np.testing.assert_array_equal(DesignSampler(pri, seed=42).sample(10), expect)
    np.testing.assert_array_equal(DesignSampler(pri, seed=42, fixed_params={}).sample(10), expect)


def test_bounds_follow_active_columns(setup):
    gl, theta, outputs, _ = setup
    emu = Emulator(glaciers=gl, explained_variance=0.999, active_params=(0, 2)).fit(theta, outputs)
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth',
                             fixed_params={'insul_scale': 1.0})
    eff = cal._effective_bounds()
    assert eff['advection_scale'] == pytest.approx((max(0.3, theta[:, 2].min()), min(1.7, theta[:, 2].max())))
    assert eff['refreeze_frac'] == pytest.approx((theta[:, 0].min(), theta[:, 0].max()))
    with pytest.raises(ValueError):   # f_ins free but not an emulator input
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth')
