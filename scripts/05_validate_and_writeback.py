#!/usr/bin/env python
"""
Step 5/5: leave-one-glacier-out validation vs Tier-2/k-NN, then write the IDL residual file.

Usage:
    python scripts/05_validate_and_writeback.py [config.yaml]

Requires scripts/04_calibrate.py to have completed (reads posterior_samples.csv via
pipeline.flat_samples, and emulator.pkl). Writes {output_dir}/heldout_points.csv,
heldout_folds.csv (the adoption test), loo_results.csv and the
firnicetemp_calibration_*_bayes_residual.dat file GloGEM's read_firnicetemp_calibration_
bayes.pro consumes -- ONLY point firnice_temp_calib_bayes_file at it if the decision rule
printed below says ADOPT and one real GloGEM run at the adopted theta confirms it.
"""

import sys

import dill
import numpy as np
import pandas as pd

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.priors import PARAM_NAMES


def main():
    config = CalibrationConfig.from_yaml(sys.argv[1]) if len(sys.argv) > 1 else CalibrationConfig()
    pipeline = CalibrationPipeline(config)
    pipeline.load_data()

    with open(config.emulator_path, 'rb') as f:
        pipeline.emulator = dill.load(f)
    if not config.posterior_path.exists():
        raise SystemExit(f'{config.posterior_path} not found -- run scripts/04_calibrate.py first.')

    # Built exactly as step 4 builds it (pipeline.make_calibrator: target, fixed parameters,
    # prior bounds, model error, points left out of the fit). Until 2026-08-31 step 5 rebuilt it
    # with defaults instead of the configured covariance, so the s_model_variance sensitivity arms
    # {2.3, 3.5, 5.0} produced BYTE-IDENTICAL loo_results.csv (md5 ef05b5d71ad5c03b40ccbd2dcf3a1984),
    # which was misread as evidence that the result is insensitive to the choice.
    pipeline.calibrator = pipeline.make_calibrator()
    if pipeline.calibrator.offset_in_fit:
        pipeline.calibrator.fit_discrepancy()
    pipeline.flat_samples = pd.read_csv(config.posterior_path)[list(PARAM_NAMES)].to_numpy()

    points, pooled, adopt = pipeline.validate_and_writeback()


if __name__ == '__main__':
    main()
