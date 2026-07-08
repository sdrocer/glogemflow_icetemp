#!/usr/bin/env python
"""
Step 4/5: fit the discrepancy GP and sample the KO posterior over theta via emcee.

Usage:
    python scripts/04_calibrate.py [config.yaml]

Requires scripts/03_fit_emulator.py to have completed (reads emulator.pkl).
Writes {output_dir}/posterior_samples.csv.
"""

import sys

import dill

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def main():
    config = CalibrationConfig.from_yaml(sys.argv[1]) if len(sys.argv) > 1 else CalibrationConfig()
    pipeline = CalibrationPipeline(config)
    pipeline.load_data()

    if not config.emulator_path.exists():
        raise SystemExit(f'{config.emulator_path} not found -- run scripts/03_fit_emulator.py first.')
    with open(config.emulator_path, 'rb') as f:
        pipeline.emulator = dill.load(f)

    pipeline.calibrate()
    print(f'posterior samples -> {config.posterior_path}')


if __name__ == '__main__':
    main()
