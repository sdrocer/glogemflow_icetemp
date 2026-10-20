"""Readings matched to their own model values by position (plan of 20 Oct, Step 2a): with one point per borehole and
year an entity can hold several points at one depth. On the synthetic setup (unique depths) every result equals the
values recorded before the change (tests/data/positional_reference.json, icetemp 8fde2b9)."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import Priors
from icetemp.calibration.validation import Validator, predict_profile_emulator

REF = json.load(open(Path(__file__).parent / 'data' / 'positional_reference.json'))
THETAS = np.array([[0.5, 1.0, 1.0], [0.3, 0.7, 1.0], [0.8, 1.4, 0.9], [0.6, 1.2, 1.1], [0.4, 0.9, 1.0]])


def test_unique_depths_unchanged(setup):
    gl, _, _, emu = setup
    cal = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth')
    np.testing.assert_allclose([cal.log_likelihood(t) for t in THETAS], REF['depth_ll'], rtol=1e-12)
    np.testing.assert_allclose([cal.mahalanobis_per_dof(t) for t in THETAS], REF['depth_maha'], rtol=1e-12)
    ent = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth')
    ent.fit_discrepancy()
    np.testing.assert_allclose([ent.log_likelihood(t) for t in THETAS], REF['entity_offset_ll'], rtol=1e-12)
    ent0 = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', offset_in_fit=False)
    np.testing.assert_allclose([ent0.log_likelihood(t) for t in THETAS], REF['entity_ll'], rtol=1e-12)
    cal2 = BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(bounds_override={'refreeze_frac': [0.1, 1.0]}),
                              track='depth', target='depth', fixed_params={'advection_scale': 1.0},
                              fit_exclude=fit_exclusions(gl, ['colle_gnifetti_below_40m']))
    points, _ = Validator(glaciers=gl, calib_df=None, calibrator=cal2).run_parameters()
    for c in ('depth', 'obs', 'pred_ko', 'pred_t1', 'pred_prior'):
        np.testing.assert_allclose(points[c], REF['points'][c], rtol=1e-12)


@pytest.fixture(scope='module')
def duplicates():
    """One entity with two holes read at the same depths, two runs with distinct values per point."""
    depths = np.array([5.0, 5.0, 12.0, 12.0, 20.0])
    g = SimpleNamespace(glacier_name='E@3005', base_glacier_name='E', glacier_id='00009', glogem_region='CentralEurope',
                        elevation=3005.0, elevation_band=3005, latitude=46.0, longitude=7.9, depths=depths,
                        T_obs=np.array([-3.0, -2.0, -2.5, -1.5, -1.0]), weights=np.ones(5), sigma=np.full(5, 0.3),
                        point_hole=np.array(['1', '2', '1', '2', '1'], dtype=object),
                        point_year=np.array([1990, 2010, 1990, 2010, 1990]))
    theta = DesignSampler(Priors(), seed=5).sample(20)
    outputs = [{g.glacier_name: (np.append(depths, 30.0),
                                 np.append(-4 + t[0] * np.arange(1, 6) + 0.3 * t[1] * np.arange(5) ** 0.5 - 0.2 * t[2],
                                           -0.5))} for t in theta]
    emu = Emulator(glaciers=[g], explained_variance=0.9999).fit(theta, outputs)
    return g, theta, outputs, emu


def test_each_point_keeps_its_model_value(duplicates):
    g, theta, outputs, emu = duplicates
    F = emu.assemble_matrix(outputs)
    for j, out in enumerate(outputs):
        np.testing.assert_array_equal(F[:5, j], out[g.glacier_name][1][:5])


def test_each_row_keeps_its_observation(duplicates):
    g, _, _, emu = duplicates
    cal = BayesianCalibrator(emulator=emu, glaciers=[g], priors=Priors(), track='depth', target='depth')
    np.testing.assert_array_equal(cal._row_obs, g.T_obs[[emu.x_points_used()[r] for r in cal._row_emu_idx]])
    assert sorted(cal._row_obs) == sorted(g.T_obs)
    val = Validator(glaciers=[g], calib_df=None, calibrator=cal)
    pred = predict_profile_emulator(emu, g.glacier_name, theta_mid := np.array([0.5, 1.0, 1.0]), g.depths,
                                    points=np.arange(5))
    mean = emu.predict(np.atleast_2d(theta_mid))[0][0]
    np.testing.assert_allclose(pred, mean[:5])
    assert len(set(np.round(pred, 9))) == 5            # no two points share a prediction


def test_excluding_a_depth_takes_every_point_there(duplicates):
    g, _, _, emu = duplicates
    cal = BayesianCalibrator(emulator=emu, glaciers=[g], priors=Priors(), track='depth', target='depth',
                             fit_exclude=frozenset({(g.glacier_name, 12.0)}))
    assert sorted(cal._row_obs) == sorted(g.T_obs[[0, 1, 4]])


def test_other_pooling_is_refused(duplicates):
    g, _, _, emu = duplicates
    other = SimpleNamespace(**{**vars(g), 'depths': np.array([5.0, 6.0, 12.0, 13.0, 20.0])})
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=[other], priors=Priors(), track='depth', target='depth')
