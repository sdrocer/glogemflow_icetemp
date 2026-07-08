#!/usr/bin/env python
"""
Step 2/5: write GloGEM training inputs, and (only with --run) launch the real IDL training
runs at every LHS design point.

Usage:
    python scripts/02_run_training.py [config.yaml]              # write inputs only (safe)
    python scripts/02_run_training.py [config.yaml] --run [idl_bin]  # actually launch IDL

Neither mode touches GloGEM/config.pro: training runs launch with the environment variable
GLOGEM_CONFIG pointed at the auto-generated training config (see runner.GloGEMRunner's module
docstring), so they are safe to run concurrently with anything else already using config.pro.
"""

import sys

import dill

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    do_run = '--run' in sys.argv
    idl_bin = args[1] if do_run and len(args) > 1 else 'idl'
    config_path = args[0] if args else None

    config = CalibrationConfig.from_yaml(config_path) if config_path else CalibrationConfig()
    pipeline = CalibrationPipeline(config)
    pipeline.load_data()
    pipeline.load_design()
    pipeline.write_training_inputs()

    if not do_run:
        print('\n(write-only mode -- pass --run [idl_bin] to actually launch IDL)')
        return

    results = pipeline.run_training(idl_bin=idl_bin)
    with open(config.output_path / 'training_results.pkl', 'wb') as f:
        dill.dump(results, f)
    print(f'training results -> {config.output_path / "training_results.pkl"}')


if __name__ == '__main__':
    main()
