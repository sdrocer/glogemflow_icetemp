"""Gather the merged five-part runs of the Step 5 design into training_results.pkl (Step 5 plan, step 2).

Design point i of <config>'s design.csv ran once per orientation part as <tag><i:02d>_FRHS_..._o<S><P>_w0
(full_set_runs.py, chain week3/chain_step5.sh) and was merged into <tag><i:02d>_FRHS_..._o<S>P_w0
(merge_part_runs.py). Checked per point before it counts: every part run holds exactly the design row's parameters
and the config's thermal settings, and the merged output covers exactly the config's calibration entities.
Writes <output_dir>/training_results.pkl ({design_NNNN: theta, ok, output}, as 02_run_training.py --merge) and
training_runs.csv (one row per point: its runs and the entities without a profile).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/assemble_training.py
           config/step5_centraleurope_parts5.yaml [--tag s5d] [--strength C] [--scheme GNESW] [--out FILE]
"""
import argparse
import pickle
import sys
from pathlib import Path

import dill
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'orientation'))
import forcing_fair_test as ft  # noqa: E402
import part_factors as pf  # noqa: E402
from icetemp.calibration import CalibrationConfig, CalibrationPipeline  # noqa: E402
from icetemp.calibration.priors import PARAM_NAMES  # noqa: E402

SETTINGS = ('forcing_toff', 'thermal_lapse', 'hl1980', 'lh_water', 'snow_refreeze', 'orient', 'thermal_air',
            'freeair_file', 'f_sur', 'catchment', 'air_temperature', 'calibration_source_dir')
# settings that runs made before them did not record, at the value those runs had
UNRECORDED = dict(air_temperature='era5_2m')


def one_run(pattern):
    hits = sorted(ft.BAYESCAL.glob(pattern))
    if len(hits) != 1:
        raise SystemExit(f'{pattern}: {len(hits)} run folders, expected 1')
    return hits[0]


def check_part(run, theta, config):
    """The part run ran the design row and the config's thermal settings."""
    d = pd.read_csv(run / 'design.csv')
    if list(d.columns) != list(PARAM_NAMES) or len(d) != 1 or not np.allclose(d.to_numpy()[0], theta, rtol=0,
                                                                              atol=1e-12):
        raise SystemExit(f'{run.name}: design {d.to_numpy().tolist()} is not {theta.tolist()}')
    c = yaml.safe_load((run / 'config.yaml').read_text())
    bad = [k for k in SETTINGS if c.get(k, UNRECORDED.get(k)) != getattr(config, k)]
    if bad:
        raise SystemExit(f'{run.name}: {", ".join(f"{k} {c.get(k)!r}" for k in bad)} differ from the config')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config')
    ap.add_argument('--tag', default='s5d')
    ap.add_argument('--strength', default='C')
    ap.add_argument('--scheme', default='GNESW', choices=sorted(pf.SCHEMES))
    ap.add_argument('--out', default=None, help='default: <output_dir>/training_results.pkl')
    a = ap.parse_args()
    config = CalibrationConfig.from_yaml(a.config)
    pl = CalibrationPipeline(config)
    pl.load_design()
    pl.load_data()
    names = {g.glacier_name for g in pl._calibration_glaciers()}
    results, rows = {}, []
    for i, theta in enumerate(pl.design):
        key = f'design_{i:04d}'
        merged = one_run(f'{a.tag}{i:02d}_*_o{a.strength}P_w0')
        parts = (merged / 'parts.txt').read_text().splitlines()[1:1 + len(pf.SCHEMES[a.scheme])]
        for p, r in zip(pf.SCHEMES[a.scheme], parts):
            if not r.endswith(f'_o{a.strength}{p}_w0'):
                raise SystemExit(f'{merged.name}/parts.txt: {r} is not part {p}')
            check_part(ft.BAYESCAL / r, theta, config)
        out = pickle.load(open(merged / 'training_runs' / 'design_0000.output_hole_year.pkl', 'rb'))
        if set(out) != names:
            raise SystemExit(f'{merged.name}: {len(set(out) ^ names)} entities differ from the config\'s '
                             f'{len(names)}')
        n_none = sum(v is None for v in out.values())
        results[key] = {'theta': theta, 'ok': True, 'output': out}
        rows.append(dict(point=key, **dict(zip(PARAM_NAMES, theta)), merged=merged.name, n_none=n_none))
    dst = Path(a.out) if a.out else config.output_path / 'training_results.pkl'
    with open(dst, 'wb') as f:
        dill.dump(results, f)
    pd.DataFrame(rows).to_csv(dst.with_name('training_runs.csv'), index=False)
    print(f'{len(results)} design points, {len(names)} entities each, '
          f'{sum(r["n_none"] for r in rows)} entity outputs without a profile -> {dst}')


if __name__ == '__main__':
    main()
