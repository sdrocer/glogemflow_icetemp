"""Validator.run_parameters / parameters_summary on the synthetic setup (conftest.py), and the
Step-5 configuration defaults."""
import numpy as np
import pytest

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors
from icetemp.calibration.validation import Validator


@pytest.fixture(scope='module')
def heldout(setup):
    gl, _, _, emu = setup
    priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=priors, track='depth', target='depth',
                             fixed_params={'advection_scale': 1.0},
                             fit_exclude=fit_exclusions(gl, ['colle_gnifetti_below_40m']))
    val = Validator(glaciers=gl, calib_df=None, calibrator=cal)
    points, folds = val.run_parameters()
    return gl, emu, cal, val, points, folds


def test_one_fold_per_system_every_point_once(heldout):
    gl, _, _, _, points, folds = heldout
    assert len(folds) == 4 and folds['n_points'].sum() == len(points)
    assert len(points) == sum(len(g.depths) for g in gl)
    assert points.groupby(['entity', 'depth']).size().max() == 1
    assert points['excluded'].sum() == 2
    assert np.all(folds['ko_advection_scale'] == 1.0) and np.all(folds['t1_advection_scale'] == 1.0)


def test_predictions_are_the_emulator_at_each_fold_theta(heldout):
    _, emu, _, _, points, folds = heldout
    x_used = emu.x_index_used()
    for _, f in folds.iterrows():
        for m in ('ko', 't1', 'prior'):
            th = [f[f'{m}_refreeze_frac'], f[f'{m}_insul_scale'], f[f'{m}_advection_scale']]
            mean = emu.predict(np.atleast_2d(th))[0][0]
            sub = points[points['system'] == f['system']]
            rows = [x_used.index((e, d)) for e, d in zip(sub['entity'], sub['depth'])]
            np.testing.assert_allclose(sub[f'pred_{m}'], mean[rows], atol=1e-12)
    assert folds['prior_refreeze_frac'].iloc[0] == pytest.approx(0.55)


def test_summary_rule(heldout):
    *_, val, points, folds = heldout
    fit = points[~points['excluded']]
    _, pooled, ok = val.parameters_summary(points, folds, margin=0.10)
    for m in ('ko', 't1', 'prior'):
        assert pooled[m] == pytest.approx(np.sqrt(np.mean((fit[f'pred_{m}'] - fit['obs']) ** 2)))
    assert ok == (pooled['ko'] < pooled['prior'] and pooled['ko'] <= pooled['t1'] + 0.10)
    _, _, ok_strict = val.parameters_summary(points, folds, margin=-10.0)
    assert not ok_strict


def test_config_defaults_round_trip(tmp_path):
    cfg = CalibrationConfig(output_dir=str(tmp_path))
    assert cfg.target == 'depth' and cfg.fixed_params == {'advection_scale': 1.0}
    assert cfg.bounds_override == {'refreeze_frac': [0.1, 1.0]} and cfg.explained_variance == 0.999
    cfg.to_yaml(tmp_path / 'c.yaml')
    back = CalibrationConfig.from_yaml(tmp_path / 'c.yaml')
    assert (back.target, back.fixed_params, back.fit_exclusions, back.bounds_override, back.rmse_margin) == \
        (cfg.target, cfg.fixed_params, cfg.fit_exclusions, cfg.bounds_override, cfg.rmse_margin)
    assert CalibrationConfig(output_dir=str(tmp_path)).fixed_params is not cfg.fixed_params
