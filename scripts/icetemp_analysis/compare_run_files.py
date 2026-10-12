"""Compare two GloGEM training runs file by file: the thermal output folder byte by byte, and the matched
model values (design_XXXX.output.pkl) exactly.

Usage:  compare_run_files.py <run_dir_a> <run_dir_b> [--pkl-a design_0000.output.pkl] [--pkl-b ...]
        [--only GLOB] [--exclude GLOB]   (e.g. --only 'firnice_snow[mb]*' for the files derived from the
        mass balance; --exclude 'firnice_wetdepth_*' against a run made before that log existed)
run_dir: data/bayescal/<tag>_w<k>. Prints identical / differing / missing files and exits 1 if any differ.
"""
import argparse
import fnmatch
import pickle
import sys
from pathlib import Path

import numpy as np

SUB = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('a')
    ap.add_argument('b')
    ap.add_argument('--pkl-a', default=None)
    ap.add_argument('--pkl-b', default=None)
    ap.add_argument('--only', default='*')
    ap.add_argument('--exclude', default=None)
    args = ap.parse_args()
    da, db = Path(args.a) / SUB, Path(args.b) / SUB
    keep = lambda n: fnmatch.fnmatch(n, args.only) and not (args.exclude and fnmatch.fnmatch(n, args.exclude))
    fa = {p.name for p in da.iterdir() if keep(p.name)}
    fb = {p.name for p in db.iterdir() if keep(p.name)}
    same = [n for n in sorted(fa & fb) if (da / n).read_bytes() == (db / n).read_bytes()]
    diff = sorted((fa & fb) - set(same))
    print(f'{len(same)} identical, {len(diff)} differ, {len(fa - fb)} only in a, {len(fb - fa)} only in b')
    for n in diff[:20]:
        print('  differs:', n)
    bad = bool(diff or fa ^ fb)
    if args.pkl_a and args.pkl_b:
        pa = pickle.load(open(Path(args.a) / 'training_runs' / args.pkl_a, 'rb'))
        pb = pickle.load(open(Path(args.b) / 'training_runs' / args.pkl_b, 'rb'))
        neq = [k for k in sorted(set(pa) | set(pb))
               if k not in pa or k not in pb or (pa[k] is None) != (pb[k] is None)
               or (pa[k] is not None and not all(np.array_equal(np.asarray(x), np.asarray(y), equal_nan=True)
                                                 for x, y in zip(pa[k], pb[k])))]
        print(f'matched values: {len(set(pa) | set(pb)) - len(neq)} entities identical, {len(neq)} differ')
        for k in neq[:10]:
            print('  differs:', k)
        bad |= bool(neq)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
