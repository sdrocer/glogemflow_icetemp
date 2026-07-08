#!/usr/bin/env python
"""
Step 5/5: leave-one-glacier-out validation vs Tier-2/k-NN, then write the IDL residual file.

Usage:
    python scripts/05_validate_and_writeback.py [config.yaml]

Requires scripts/04_calibrate.py to have completed (reads posterior_samples.csv via
pipeline.flat_samples, and emulator.pkl). Writes {output_dir}/loo_results.csv and the
firnicetemp_calibration_*_bayes_residual.dat file GloGEM's read_firnicetemp_calibration_
bayes.pro consumes -- ONLY point firnice_temp_calib_bayes_file at it if the decision rule
printed below says ADOPT.
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

    from icetemp.calibration.calibrator import BayesianCalibrator
    from icetemp.calibration.priors import Priors
    pipeline.priors = Priors()
    pipeline.calibrator = BayesianCalibrator(
        emulator=pipeline.emulator, glaciers=pipeline._calibration_glaciers(), priors=pipeline.priors,
    )
    pipeline.calibrator.fit_discrepancy()
    pipeline.flat_samples = pd.read_csv(config.posterior_path)[list(PARAM_NAMES)].to_numpy()

    results, means, adopt = pipeline.validate_and_writeback()


if __name__ == '__main__':
    main()
