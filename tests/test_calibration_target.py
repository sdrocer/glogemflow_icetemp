"""BayesianCalibrator targets on the synthetic setup (conftest.py): the depth target's Woodbury
likelihood and covariance, the entity target's covariance with and without the offset, and points
left out of the fit."""
import numpy as np
import pytest

from icetemp.calibration.calibrator import BayesianCalibrator, _gaussian_logpdf
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors

THETAS = np.array([[0.5, 1.0, 1.0], [0.3, 0.7, 1.0], [0.8, 1.4, 0.9], [0.6, 1.2, 1.1], [0.4, 0.9, 1.0]])


def test_depth_target_fast_equals_full(setup):
    gl, _, _, emu = setup
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth')
    for th in THETAS:
        fast = cal.log_likelihood(th)
        full = _gaussian_logpdf(*cal.residual_and_covariance(th))
        assert np.isfinite(fast)
        assert abs(fast - full) <= 1e-9 * abs(full)
        r, S = cal.residual_and_covariance(th)
        assert cal.mahalanobis_per_dof(th) == pytest.approx(r @ np.linalg.solve(S, r) / len(r), rel=1e-9)


def test_depth_target_rows_and_covariance(setup):
    gl, _, _, emu = setup
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth')
    n = sum(len(g.depths) for g in gl)
    th = THETAS[0]
    r, S = cal.residual_and_covariance(th)
    assert len(r) == n and S.shape == (n, n)
    mean, L = emu.predict_factor(th)
    rows = cal._row_emu_idx
    obs = np.concatenate([g.T_obs for g in gl])
    np.testing.assert_allclose(r, obs - mean[rows], rtol=0, atol=1e-12)
    z = np.concatenate([g.depths for g in gl])
    ent = np.repeat(np.arange(len(gl)), [len(g.depths) for g in gl])
    S_model = np.where(ent[:, None] == ent[None, :], 3.5 * np.exp(-np.abs(z[:, None] - z[None, :]) / 15.0), 0)
    expect = np.diag(0.3 ** 2 + emu.truncation_variance() + np.zeros(n)) + S_model + L[rows] @ L[rows].T
    np.testing.assert_allclose(S, expect, rtol=1e-12, atol=1e-12)


def test_entity_target_covariance(setup):
    gl, _, _, emu = setup
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth')
    assert cal.target == 'entity' and cal.offset_in_fit
    cal.fit_discrepancy()
    th = THETAS[2]
    r, S = cal.residual_and_covariance(th)
    mean, L = emu.predict_factor(th)
    A, rows = cal._A, cal._row_emu_idx
    M = np.diag(cal._row_sigma_obs2 + emu.truncation_variance()) + cal._row_model_covariance()
    np.testing.assert_allclose(r, A @ (cal._row_obs - mean[rows]), atol=1e-12)
    np.testing.assert_allclose(S, cal._K_delta + A @ M @ A.T + (A @ L[rows]) @ (A @ L[rows]).T, atol=1e-12)
    no_offset = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', offset_in_fit=False)
    r2, S2 = no_offset.residual_and_covariance(th)   # no fit_discrepancy needed
    np.testing.assert_allclose(S2, S - cal._K_delta, atol=1e-12)
    np.testing.assert_array_equal(r2, r)


def test_depth_target_refuses_offset_and_basal(setup):
    gl, _, _, emu = setup
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth',
                           offset_in_fit=True)
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='all', target='depth')


def test_fit_exclusions(setup):
    gl, _, _, emu = setup
    excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
    assert excl == {('C@4455', 45.0), ('C@4455', 55.0)}
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth',
                             fit_exclude=excl)
    assert len(cal._row_obs) == sum(len(g.depths) for g in gl) - 2
    assert fit_exclusions(gl, []) == frozenset()
    with pytest.raises(ValueError):
        fit_exclusions(gl, ['no_such_rule'])
