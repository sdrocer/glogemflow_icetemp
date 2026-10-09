"""Re-match an ensemble's archived raw profiles with the current comparison, without IDL.

Usage:  rematch_ensemble.py <tag> <out_dir>
Reads net_scratch/icetemp_profiles/<tag>_w<k>/design_XXXX.tar.gz (k = 0..5) for the entities in the
stored outputs of data/bayescal/<tag>_w<k>/training_runs and writes <out_dir>/design_XXXX.pkl
(dict entity -> (depths, temperatures)), the format score_paired_ensembles.py reads.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/rematch_ensemble.py
"""
import pickle
import sys
from pathlib import Path

from icetemp.calibration.data import DataHandler
from icetemp.calibration.runner import GloGEMRunner

REPO = Path(__file__).resolve().parents[2]


def main():
    tag, out = sys.argv[1], Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    dh = DataHandler(region='CentralEurope', include_estimated=True, include_unassessed=True)
    dh.load()
    by_name = {g.glacier_name: g for g in dh.calibration_glaciers}
    for w in range(6):
        run_dir = REPO / 'data' / 'bayescal' / f'{tag}_w{w}' / 'training_runs'
        runner = GloGEMRunner(run_dir=run_dir)
        for f in sorted(run_dir.glob('design_*.output.pkl')):
            name = f.name.split('.')[0]
            stored = pickle.load(open(f, 'rb'))
            gs = [by_name[n] for n in stored if n in by_name]
            if not runner.profile_archive(name).exists():
                print(f'{name}: no archive, skipped', flush=True)
                continue
            pickle.dump(runner.reparse_archived(gs, name), open(out / f'{name}.pkl', 'wb'))
            print(f'{name}: {len(gs)} entities', flush=True)


if __name__ == '__main__':
    main()
