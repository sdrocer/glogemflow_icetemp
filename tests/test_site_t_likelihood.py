"""Form C6b of the rules plan (27 Oct): the per-reading term at shallow depth, the multivariate t per site
against scipy, its Gaussian limit, site weights, the guards, and the sizes recovered from simulated data."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import multivariate_normal, multivariate_t

from icetemp.calibration import error_model as em
from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import Priors
from icetemp.calibration.residual_covariance import model_covariance

THETAS = np.array([[0.5, 1.0, 1.0], [0.3, 0.7, 1.0], [0.8, 1.4, 0.9]])
SIZES = dict(s_model_variance=2.3, s_model_length_m=60.0, s_hole_variance=1.4, s_hole_year_variance=0.05,
             s_site_variance=0.4, s_white_variance=0.37)


def entity(name, holes, years, depths, rng):
    depths = np.asarray(depths, dtype=float)
    return SimpleNamespace(glacier_name=name, base_glacier_name=name.split('@')[0], glacier_id='0000' + name[0],
                           glogem_region='CentralEurope', elevation=float(name.split('@')[1]),
                           elevation_band=int(name.split('@')[1]), latitude=46.0, longitude=7.9, depths=depths,
                           T_obs=-3 + 0.05 * depths + rng.normal(0, 1.5, len(depths)), weights=np.ones(len(depths)),
                           sigma=np.full(len(depths), 0.2), point_hole=np.array(holes, dtype=object),
                           point_year=np.array(years))


@pytest.fixture(scope='module')
def three_sites():
    """Sites P (two bands, not neighbours in the list), Q and R; depths on both sides of 15 m."""
    rng = np.random.default_rng(5)
    gl = [entity('P@3005', ['1', '1', '2', '1', '2'], [1990, 2010, 2010, 1990, 2010], [5, 5, 5, 22, 22], rng),
          entity('Q@2505', ['7', '7', '7', '7'], [2000, 2000, 2000, 2000], [4, 9, 20, 35], rng),
          entity('P@3015', ['3', '3', '3'], [2005, 2005, 2005], [3, 8, 16], rng),
          entity('R@3505', ['9', '9', '8'], [2001, 2001, 2001], [2, 18, 30], rng)]
    theta = DesignSampler(Priors(), seed=4).sample(20)
    outputs = [{g.glacier_name: (g.depths, -5 + t[0] * np.exp(-g.depths / 10) + 0.4 * t[1] - 0.2 * t[2]
                                 + 0.01 * g.depths) for g in gl} for t in theta]
    return gl, Emulator(glaciers=gl, explained_variance=0.9999).fit(theta, outputs)


def calibrator(gl, emu, **kw):
    return BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), track='depth', target='depth',
                              site_of={g.glacier_name: g.glacier_name.split('@')[0] for g in gl}, **kw)


def per_site_reference(cal, theta, nu=None):
    """scipy's log density of each site's rows under the dense covariance's site block."""
    r, S = cal.residual_and_covariance(theta)
    site = np.asarray(cal._row_site, dtype=object)
    out = []
    for s in cal._block_site:
        i = np.flatnonzero(site == s)
        Sk = S[np.ix_(i, i)]
        if nu is None:
            out.append(multivariate_normal(mean=np.zeros(len(i)), cov=Sk).logpdf(r[i]))
        else:
            out.append(multivariate_t(loc=np.zeros(len(i)), shape=Sk * (nu - 2) / nu, df=nu).logpdf(r[i]))
    return np.array(out)


def test_white_term_only_at_shallow_depth():
    C = model_covariance(['E', 'E', 'E'], [5.0, 15.0, 15.5], s2=0.0, ell=1.0, w2=0.4, w_max_depth=15.0)
    np.testing.assert_allclose(C, np.diag([0.4, 0.4, 0.0]))
    np.testing.assert_allclose(model_covariance(['E', 'F'], [5.0, 40.0], s2=0.0, ell=1.0, w2=0.4), 0.4 * np.eye(2))


@pytest.mark.parametrize('error', ['mean', 'varying'])
def test_t_per_site_equals_scipy(three_sites, error):
    gl, emu = three_sites
    cal = calibrator(gl, emu, site_t_nu=2.6, emulator_error=error, **SIZES)
    for th in THETAS:
        ref = per_site_reference(cal, th, nu=2.6)
        np.testing.assert_allclose(cal.site_log_likelihoods(th), ref, rtol=1e-10)
        assert abs(cal.log_likelihood(th) - ref.sum()) <= 1e-9 * abs(ref.sum())


def test_t_tends_to_the_gaussian(three_sites):
    gl, emu = three_sites
    gauss = calibrator(gl, emu, **SIZES)
    t_big = calibrator(gl, emu, site_t_nu=1e8, **SIZES)
    for th in THETAS:
        assert abs(t_big.log_likelihood(th) - gauss.log_likelihood(th)) < 1e-4


def test_site_weights(three_sites):
    gl, emu = three_sites
    cal = calibrator(gl, emu, site_weights={'P': 0.5}, **SIZES)
    w = np.array([0.5 if s == 'P' else 1.0 for s in cal._block_site])
    for th in THETAS:
        ref = per_site_reference(cal, th)
        assert abs(cal.log_likelihood(th) - w @ ref) <= 1e-9 * abs(w @ ref)
    with pytest.raises(ValueError):
        calibrator(gl, emu, site_weights={'nowhere': 0.5}, **SIZES)


def test_guards(three_sites):
    gl, emu = three_sites
    with pytest.raises(ValueError):
        calibrator(gl, emu, site_t_nu=2.0, **SIZES)
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), target='entity', s_white_variance=0.3)
    with pytest.raises(ValueError):
        BayesianCalibrator(emulator=emu, glaciers=gl, priors=Priors(), target='entity', site_t_nu=3.0)


def test_estimator_likelihood_equals_scipy(three_sites):
    gl, emu = three_sites
    cal = calibrator(gl, emu, site_t_nu=3.0, emulator_error='mean', **SIZES)
    th = THETAS[0]
    p = em.sizes_of(cal)
    blocks = em.SiteBlocks(cal)
    r = em.residual(cal, th)
    assert abs(-em.neg_log_likelihood(blocks, r, p) - per_site_reference(cal, th, nu=3.0).sum()) < 1e-8
    assert abs(-em.neg_log_likelihood(blocks, r, p, gaussian=True) - per_site_reference(cal, th).sum()) < 1e-8
    np.testing.assert_allclose(r, cal.residual_and_covariance(th)[0])


def simulated_structure(n_sites=60, seed=0):
    """Rows of a fake calibrator: per site two bands, per band two holes (one read in two years), 12 depths."""
    rng = np.random.default_rng(seed)
    ent, dep, hole, year, site = [], [], [], [], []
    depths = np.array([1, 3, 5, 8, 11, 14, 17, 22, 28, 35, 45, 60], float)
    for s in range(n_sites):
        for b in range(2):
            for h in range(2):
                for y in ((1990, 2010) if h == 0 else (2000,)):
                    k = len(set(ent))
                    ent += [f'{s}@{b}'] * len(depths)
                    dep += list(depths + rng.uniform(-0.4, 0.4, len(depths)))
                    hole += [f'{s}-{b}-{h}'] * len(depths)
                    year += [y] * len(depths)
                    site += [f'site{s}'] * len(depths)
    names = {e: k for k, e in enumerate(dict.fromkeys(ent))}
    return SimpleNamespace(target='depth', _rows_ok=True, emulator=SimpleNamespace(truncation_variance=lambda: 0.0),
                           s_white_max_depth=15.0, _row_site=site, _row_entity=np.array([names[e] for e in ent]),
                           _row_depth=np.array(dep), _row_hole=hole, _row_year=year,
                           _row_sigma_obs2=np.full(len(dep), 0.04))


def test_sizes_recovered_from_simulated_data():
    cal = simulated_structure()
    truth = dict(s2=2.0, ell=50.0, sh2=1.0, sy2=0.3, ss2=0.5, w2=0.4, nu=4.0)
    blocks = em.SiteBlocks(cal)
    rng = np.random.default_rng(11)
    r = np.empty(blocks.n)
    c = (truth['nu'] - 2) / truth['nu']
    for i in blocks.idx:
        z = rng.multivariate_normal(np.zeros(len(i)), blocks.covariance(truth, i))
        r[i] = z * np.sqrt(c * truth['nu'] / rng.chisquare(truth['nu']))
    fit = em.fit_sizes(cal, None, r=r)
    got = {n: fit[em.FIELD[n]] for n in em.NAMES}
    for n, tol in dict(s2=0.35, ell=0.6, sh2=0.5, ss2=0.8, w2=0.25).items():
        assert abs(np.log(got[n] / truth[n])) < tol, (n, got[n], truth[n])
    assert 2.5 < got['nu'] < 12, got['nu']
    assert -fit['nll'] >= -em.neg_log_likelihood(blocks, r, truth) - 1e-6
