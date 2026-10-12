#!/usr/bin/env python
"""
Step 1b: the second design stage. Appends config.stage2_points runs around the stage-1 result
(grid posterior mean +- stage2_half_width_sd posterior sd, written by scripts/04_calibrate.py) to
design.csv.

Usage:
    python scripts/01b_add_stage2_design.py [config.yaml]

Then scripts/02_run_training.py --run runs only the new rows, and steps 3-5 refit on all of them.
"""

import sys

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def main():
    config = CalibrationConfig.from_yaml(sys.argv[1]) if len(sys.argv) > 1 else CalibrationConfig()
    CalibrationPipeline(config).add_stage2_design()


if __name__ == '__main__':
    main()
