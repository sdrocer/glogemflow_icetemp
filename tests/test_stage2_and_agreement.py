"""The second design stage around a first result, and the agreement check of the fits with and
without the emulator's uncertainty."""
import json

import numpy as np
import pandas as pd

from icetemp.calibration import adoption
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.pipeline import CalibrationPipeline
from icetemp.calibration.priors import PARAM_NAMES, Priors


def test_sample_around_stays_in_its_box():
    pri = Priors(bounds_override={'refreeze_frac': (0.1, 1.0)})
    s = DesignSampler(pri, seed=5, fixed_params={'advection_scale': 1.0})
    d = s.sample_around(20, np.array([0.15, 0.8, 9.0]), np.array([0.1, 0.25, 9.0]))
    assert d.shape == (20, 3) and np.all(d[:, 2] == 1.0)
    assert d[:, 0].min() >= 0.1 and d[:, 0].max() <= 0.25      # clipped at the prior bound 0.1
    assert d[:, 1].min() >= 0.55 and d[:, 1].max() <= 1.05
    # one point per stratum: a Latin hypercube
    assert len(np.unique(np.floor((d[:, 1] - 0.55) / 0.5 * 20))) == 20
    np.testing.assert_array_equal(d, s.sample_around(20, np.array([0.15, 0.8, 9.0]), np.array([0.1, 0.25, 9.0])))


def test_add_stage2_design_appends(tmp_path):
    cfg = CalibrationConfig(output_dir=str(tmp_path), n_design_points=40, stage2_points=20)
    pl = CalibrationPipeline(cfg)
    pl.priors = Priors(bounds_override=cfg.bounds_override)
    first = DesignSampler(pl.priors, seed=cfg.design_seed, fixed_params=cfg.fixed_params).sample_df(40)
    first.to_csv(cfg.design_path, index=False)
    before = cfg.design_path.read_text()
    cfg.diagnostics_path.write_text(json.dumps(dict(grid_names=['refreeze_frac', 'insul_scale'],
                                                    grid_mean=[0.35, 0.78], grid_sd=[0.03, 0.12])))
    new = pl.add_stage2_design()
    assert cfg.design_path.read_text().startswith(before)       # stage-1 rows untouched
    d = pd.read_csv(cfg.design_path, float_precision='round_trip')[list(PARAM_NAMES)].to_numpy()
    assert d.shape == (60, 3)
    np.testing.assert_array_equal(d[:40], first.to_numpy())
    np.testing.assert_array_equal(d[40:], new)
    assert np.all(np.abs(new[:, 0] - 0.35) <= 2.5 * 0.03 + 1e-12)
    assert np.all(np.abs(new[:, 1] - 0.78) <= 2.5 * 0.12 + 1e-12) and np.all(new[:, 2] == 1.0)


def moments(mean, sd, corr):
    return dict(names=('refreeze_frac', 'insul_scale'), mean=np.array(mean), sd=np.array(sd),
                corr=np.array([[1.0, corr], [corr, 1.0]]))


def test_agreement_check():
    prior_sd = {'refreeze_frac': 0.26, 'insul_scale': 0.33}
    base = moments([0.35, 0.78], [0.03, 0.12], 0.65)
    ok = adoption.check_emulator_agreement(base, moments([0.36, 0.80], [0.035, 0.13], 0.7), prior_sd)
    assert ok.passed and not ok.reasons
    # the Step 5b trial on alps_s2: with E the fit blurs (grid: 0.32/0.62, sd 0.086/0.287, corr 0.96)
    trial = adoption.check_emulator_agreement(base, moments([0.319, 0.621], [0.086, 0.287], 0.96), prior_sd)
    assert not trial.passed and any('correlate' in r for r in trial.reasons)
    assert any('insul_scale keeps' in r for r in trial.reasons)
    moved = adoption.check_emulator_agreement(base, moments([0.45, 0.78], [0.03, 0.12], 0.6), prior_sd)
    assert not moved.passed and moved.shift_sd['refreeze_frac'] > 1.0
