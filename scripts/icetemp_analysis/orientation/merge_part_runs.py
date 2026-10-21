"""Merge the three orientation part runs of one parameter set into one run (Step 5 plan, step 1c).

The runs <base>_o<S>G_w0, <base>_o<S>N_w0, <base>_o<S>S_w0 (full_set_runs.py --orient-dir .../orient_parts/<S>_<P>
--orient-tag <S><P>) are rematched per borehole and year from their archived profiles (as rescore_hole_year.py),
merged by runner.GloGEMRunner.merge_part_outputs with the hole parts of hole_parts.csv, and written as the run
<base>_o<S>P_w0: training_runs/design_0000.output_hole_year.pkl, plus parts.txt naming the three source runs, so the
scorers read the merged run like any other (pooled per hole and year; the merge is defined per hole, so there is no
entity-pooled output; profiles and frozen beds are read from the part run of each hole). The mass-balance files of
the three runs must be identical (snow logs and firn line), which is checked first.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/orientation/merge_part_runs.py
           BASE [BASE ...] [--strength C] [--parts DIR]
       BASE e.g. o3_FRHS_r0.85_i1.017_a1.00
"""
import argparse
import filecmp
import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration.data import DataHandler  # noqa: E402
from icetemp.calibration.runner import GloGEMRunner  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
PARTS = ('G', 'N', 'S')


def same_mass_balance(runs):
    """The files derived from the mass balance are identical in the three part runs."""
    d = [ft.BAYESCAL / r / ft.FIRNICE for r in runs]
    names = sorted(p.name for p in d[0].glob('firnice_snow[bm]*_*.dat')) + \
        sorted(p.name for p in d[0].glob('firnice_firnline_*.dat'))
    if not names:
        raise SystemExit(f'{runs[0]}: no mass-balance files to compare')
    bad = [n for n in names for x in d[1:] if not filecmp.cmp(d[0] / n, x / n, shallow=False)]
    if bad:
        raise SystemExit(f'mass-balance files differ between the part runs: {bad[:5]}')
    return len(names)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('bases', nargs='+')
    ap.add_argument('--strength', default='C')
    ap.add_argument('--parts', default=str(W3 / 'orient_parts'))
    a = ap.parse_args()
    hp = pd.read_csv(Path(a.parts) / 'hole_parts.csv', dtype={'hole': str})
    part_of_hole = dict(zip(hp.hole, hp.part))
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    for base in a.bases:
        runs = [f'{base}_o{a.strength}{p}_w0' for p in PARTS]
        n_mb = same_mass_balance(runs)
        outs = {}
        for p, r in zip(PARTS, runs):
            runner = GloGEMRunner(run_dir=ft.BAYESCAL / r / 'training_runs')
            stored = pickle.load(open(ft.output_path(r), 'rb'))
            gl = [g for g in dh.calibration_glaciers if g.glacier_name in stored]
            outs[p] = runner.reparse_archived(gl, 'design_0000')
        merged = GloGEMRunner.merge_part_outputs(outs, gl, part_of_hole)
        dst = ft.BAYESCAL / f'{base}_o{a.strength}P_w0' / 'training_runs'
        dst.mkdir(parents=True, exist_ok=True)
        with open(dst / 'design_0000.output_hole_year.pkl', 'wb') as f:
            pickle.dump(merged, f)
        (dst.parent / 'parts.txt').write_text('merged from the part runs\n' + '\n'.join(runs) + '\n'
                                              f'hole parts: {Path(a.parts) / "hole_parts.csv"}\n')
        n_none = sum(v is None for v in merged.values())
        print(f'{base}: merged {len(merged)} entities ({n_none} None), {n_mb} mass-balance files identical',
              flush=True)


if __name__ == '__main__':
    main()
