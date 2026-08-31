#!/usr/bin/env python
"""
Step 4/5: fit the discrepancy GP and sample the KO posterior over theta via emcee.

Usage:
    python scripts/04_calibrate.py [config.yaml]

Requires scripts/03_fit_emulator.py to have completed (reads emulator.pkl).
Writes {output_dir}/posterior_samples.csv.
"""

import os

# PIN BLAS THREADS BEFORE NUMPY IS IMPORTED -- must precede any numpy/scipy import to take
# effect. The likelihood's matrices are small (107x107, or 695x695 for the depth-resolved
# variant); BLAS threading buys nothing on them and actively hurts on a shared box. Measured
# 2026-08-20 with 22 concurrent GMIP4 IDL runs: pinned gives 5.6-8.9 ms/eval with a tight tail,
# while UNPINNED was the only setting that produced 100 ms - 2188 ms outliers from
# oversubscription -- and is almost certainly the source of a spurious "1.1 s/eval, 293 h of
# MCMC" figure that briefly looked like a blocking objection. 4 measured fastest and tightest;
# 1 is the good-citizen choice alongside other users' jobs.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '4')

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
