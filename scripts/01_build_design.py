#!/usr/bin/env python
"""
Step 1/5 of the Tier-3 Bayesian calibration: load glenglat data, build the LHS design.

Usage:
    python scripts/01_build_design.py [config.yaml]

Writes {output_dir}/design.csv (the theta = (perm_frac, dT_scale, z0) design matrix). Does
NOT touch IDL or config.pro.
"""

import sys

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def main():
    config = CalibrationConfig.from_yaml(sys.argv[1]) if len(sys.argv) > 1 else CalibrationConfig()
    config.to_yaml(config.output_path / 'config.yaml')

    pipeline = CalibrationPipeline(config)
    pipeline.load_data()
    pipeline.build_design()


if __name__ == '__main__':
    main()
