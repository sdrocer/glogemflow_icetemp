#!/usr/bin/env python
"""
Step 3/5: fit the SVD + surmise PCGPwM emulator on completed training runs.

Usage:
    python scripts/03_fit_emulator.py [config.yaml]

Requires scripts/02_run_training.py --run to have completed (reads training_results.pkl).
Writes {output_dir}/emulator.pkl (dill-serialized -- plain pickle cannot serialize surmise's
internal emulator object; dill is already a surmise dependency, see requirements check in
this project's session notes).
"""

import sys

import dill

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def main():
    config = CalibrationConfig.from_yaml(sys.argv[1]) if len(sys.argv) > 1 else CalibrationConfig()
    pipeline = CalibrationPipeline(config)
    pipeline.load_data()
    pipeline.load_design()

    results_path = config.output_path / 'training_results.pkl'
    if not results_path.exists():
        raise SystemExit(
            f'{results_path} not found -- run scripts/02_run_training.py --run first.'
        )
    with open(results_path, 'rb') as f:
        training_results = dill.load(f)

    pipeline.fit_emulator(training_results)
    with open(config.emulator_path, 'wb') as f:
        dill.dump(pipeline.emulator, f)
    print(f'emulator -> {config.emulator_path}')


if __name__ == '__main__':
    main()
