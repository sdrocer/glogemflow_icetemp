"""The borehole, borehole-year and site terms of the depth target (plan of 20 Oct, Step 2b): the covariance
formula, the site-blocked factorisation against the dense likelihood, the guards, and the site table."""
from types import SimpleNamespace

import numpy as np
import pytest

from icetemp.calibration.calibrator import BayesianCalibrator, _gaussian_logpdf
from icetemp.calibration.data import site_keys
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import Priors
from icetemp.calibration.residual_covariance import model_covariance

THETAS = np.array([[0.5, 1.0, 1.0], [0.3, 0.7, 1.0], [0.8, 1.4, 0.9]])


def test_formula_on_the_example():
    """The three readings at 16 m of chapter 11: hole A 1991, hole A 2003, hole B 2018."""
    C = model_covariance(['E', 'E', 'E'], [16.0, 16.0, 16.0], s2=3.0, ell=15.0, hole=['A', 'A', 'B'],
                         year=[1991, 2003, 2018], site=['S', 'S', 'S'], sh2=2.0, sy2=0.5, ss2=1.0)
    expect = 3.0 + 1.0 + np.array([[2.0 + 0.5, 2.0, 0.0], [2.0, 2.0 + 0.5, 0.0], [0.0, 0.0, 2.0 + 0.5]])
    np.testing.assert_allclose(C, expect)
    np.testing.assert_array_equal(model_covariance(['E', 'F'], [1.0, 2.0], s2=3.0, ell=15.0),
                                  np.where(np.eye(2) > 0, 3.0, 0.0))
    with pytest.raises(ValueError):
        model_covariance(['E', 'E'], [1.0, 2.0], s2=1.0, ell=15.0, hole=[None, None], year=[None, None], sh2=1.0)


def entity(name, holes, years, depths, rng):
    depths = np.asarray(depths, dtype=float)
    return SimpleNamespace(glacier_name=name, base_glacier_name=name.split('@')[0], glacier_id='0000' + name[0],
                           glogem_region='CentralEurope', elevation=float(name.split('@')[1]),
                           elevation_band=int(name.split('@')[1]), latitude=46.0, longitude=7.9, depths=depths,
                           T_obs=-3 + 0.05 * depths + rng.normal(0, 0.5, len(depths)), weights=np.ones(len(depths)),
                           sigma=np.full(len(depths), 0.2), point_hole=np.array(holes, dtype=object),
                           point_year=np.array(years))


@pytest.fixture(scope='module')
def two_sites():
    """Site P with two bands (P@3005 two holes, one read twice; P@3015), site Q with one band; the
    entities of P are not neighbours in the list, so the site blocks need the permutation."""
    rng = np.random.default_rng(2)
    gl = [entity('P@3005', ['1', '1', '2', '1', '2'], [1990, 2010, 2010, 1990, 2010], [5, 5, 5, 12, 12], rng),
          entity('Q@2505', ['7', '7', '7'], [2000, 2000, 2000], [4, 9, 20], rng),
          entity('P@3015', ['3', '3', '3'], [2005, 2005, 2005], [3, 8, 16], rng)]
    theta = DesignSampler(Priors(), seed=4).sample(20)
    outputs = [{g.glacier_name: (g.depths, -5 + t[0] * np.exp(-g.depths / 10) + 0.4 * t[1] - 0.2 * t[2]
                                 + 0.01 * g.depths) for g in gl} for t in theta]
    return gl, Emulator(glaciers=gl, explained_variance=0.9999).fit(theta, outputs)


def calibrator(gl, emu, **kw):
    return BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth',
                              site_of={g.glacier_name: g.glacier_name.split('@')[0] for g in gl}, **kw)


@pytest.mark.parametrize('error', ['mean', 'varying'])
def test_site_blocks_equal_the_dense_likelihood(two_sites, error):
    gl, emu = two_sites
    cal = calibrator(gl, emu, s_hole_variance=2.0, s_hole_year_variance=0.5, s_site_variance=1.0,
                     emulator_error=error)
    for th in THETAS:
        fast = cal.log_likelihood(th)
        full = _gaussian_logpdf(*cal.residual_and_covariance(th))
        assert np.isfinite(fast) and abs(fast - full) <= 1e-9 * abs(full)


def test_covariance_holds_every_term(two_sites):
    gl, emu = two_sites
    cal = calibrator(gl, emu, s_hole_variance=2.0, s_hole_year_variance=0.5, s_site_variance=1.0)
    M = cal._row_model_covariance()
    e, h, y, s = cal._row_entity, np.array(cal._row_hole), np.array(cal._row_year), np.array(cal._row_site)
    d = cal._row_depth
    expect = (np.where(e[:, None] == e[None, :], 3.5 * np.exp(-np.abs(d[:, None] - d[None, :]) / 15.0), 0)
              + 2.0 * (h[:, None] == h[None, :]) + 0.5 * ((h[:, None] == h[None, :]) & (y[:, None] == y[None, :]))
              + 1.0 * (s[:, None] == s[None, :]))
    np.testing.assert_allclose(M, expect, atol=1e-12)
    zero = calibrator(gl, emu)
    no_terms = model_covariance(zero._row_entity, zero._row_depth, 3.5, 15.0)
    np.testing.assert_array_equal(zero._row_model_covariance(), no_terms)


def test_guards(two_sites):
    gl, emu = two_sites
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', s_hole_variance=1.0)
    pooled = [SimpleNamespace(**{**vars(g), 'point_hole': None, 'point_year': None}) for g in gl]
    with pytest.raises(ValueError):
        calibrator(pooled, emu, s_hole_variance=1.0)


def test_site_table():
    gl = [SimpleNamespace(glacier_name=n, elevation_band=int(n.split('@')[1].split('~')[0])) for n in
          ('Grenzgletscher@2605', 'Grenzgletscher@4455~01225', 'Taconnaz Glacier@3415', 'Glacier de Taconnaz@3835',
           'Glacier de Tête Rousse@3135')]
    s = site_keys(gl)
    assert s['Grenzgletscher@2605'] != s['Grenzgletscher@4455~01225']
    assert s['Taconnaz Glacier@3415'] == s['Glacier de Taconnaz@3835']
    assert s['Glacier de Tête Rousse@3135'] == 'Glacier de Tête Rousse'
