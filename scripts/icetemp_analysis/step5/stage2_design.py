"""Step 5, step 4.1: the 20 further design points where the answer lies (plan of 20 Oct; width rule fixed 28 Oct).

Box: centre the main fit's result (fit_step5.py posterior mean; f_ins 1.0 if rule 4 fixed it); half-width per
parameter of the joint fit the larger of HALF_SD times the result's own uncertainty (step5_common.result_sd) and the
largest distance of a fold's result (benchmark.py folds) from the centre, the same on both sides, so that the box holds
where the left-out massifs put the answer; clipped to the prior bounds. f_ins varies in the box also when rule 4 fixed
it, so that the refit on all 60 runs decides rule 4 again. Points by DesignSampler.sample_around (seed design_seed + 1,
the config's fixed parameters), appended to design.csv as rows 40-59 (run tags s5d40-s5d59, the chain of the first 40).
Refuses unless the fit and the folds were computed on the current training set of the first 40 runs.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/stage2_design.py [CONFIG]
           [--fit DIR] [--folds DIR] [--dry] [--test]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step5_common as sc  # noqa: E402
from icetemp.calibration import CalibrationConfig  # noqa: E402
from icetemp.calibration.design import DesignSampler  # noqa: E402
from icetemp.calibration.priors import PARAM_NAMES, Priors  # noqa: E402

HALF_SD = 2.5


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--fit', default=None)
    ap.add_argument('--folds', default=None)
    ap.add_argument('--dry', action='store_true', help='print the box and points, write nothing')
    ap.add_argument('--test', action='store_true', help='a code test: skips the committed-tree check')
    a = ap.parse_args()
    head = sc.rules_commit(a.config, test=a.test)
    config = CalibrationConfig.from_yaml(a.config)
    F = json.loads((Path(a.fit or config.output_path / 'fit') / 'fit.json').read_text())
    B = json.loads((Path(a.folds or config.output_path / 'benchmark') / 'folds.json').read_text())
    design = pd.read_csv(config.design_path)
    if len(design) != config.n_design_points:
        raise SystemExit(f'design.csv has {len(design)} rows, expected the first {config.n_design_points} only')
    train = dict(n_runs=config.n_design_points, sha256=sc.sha256(config.output_path / 'training_results.pkl'))
    for what, rec in (('fit.json', F['training']), ('folds.json', B['training'])):
        if rec != train:
            raise SystemExit(f'{what} was computed on {rec}, not on the current training set {train}')
    folds = B['folds']
    result = np.asarray(F['posterior_mean'])
    sd = sc.result_sd(F)
    centre, half = result.copy(), np.zeros(len(PARAM_NAMES))
    for name in F['grid']['names']:
        k = PARAM_NAMES.index(name)
        reach = max(abs(f['theta'][k] - result[k]) for f in folds)
        half[k] = max(HALF_SD * sd[name], reach)
    priors = Priors(bounds_override=config.bounds_override)
    new = DesignSampler(priors=priors, seed=config.design_seed + 1,
                        fixed_params=F['config_fixed_params']).sample_around(config.stage2_points, centre, half)
    print(f'rules commit {head}')
    print(f'centre {np.round(centre, 4).tolist()} half-width {np.round(half, 4).tolist()}')
    print(pd.DataFrame(new, columns=PARAM_NAMES).round(4).to_string())
    if not a.dry:
        pd.DataFrame(new, columns=PARAM_NAMES).to_csv(config.design_path, mode='a', header=False, index=False)
        print(f'appended {len(new)} rows -> {config.design_path}')


if __name__ == '__main__':
    main()
