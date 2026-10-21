"""Burial tables for the perfect-burial test (GloGEM firnice_burial_dir; burial_rule.py): the measured long-term net
accumulation at the borehole bands that have one, keyed by the model band GloGEM writes the hole's profile at
(band_model of the orientation investigation's band_forcing.csv). Where two entities share a model band, the band
takes their measured values weighted by scored points (as write_orient_tables.py). Glaciers without a measured band
get no file and keep the model's b_ref.

Input: a CSV with columns entity, acc_central [m w.e./yr], confidence, reference (the literature compilation of
21 Oct, week3/burial_research/accumulation.csv). Writes OUT/M/burial_<id>.dat and OUT/burial.csv (per entity: the
model's b_ref from the bref run's log, the measured value, the value applied to its band).

Usage: PYTHONPATH=src python scripts/icetemp_analysis/burial/write_burial_tables.py ACC_CSV [OUT]
"""
import re
import sys
from pathlib import Path

import pandas as pd

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
ORI = W3 / 'orientation_investigation'
BREF_LOG = (Path(__file__).resolve().parents[3] / 'data/bayescal/bref_FRHS_r0.85_i1.017_a1.00_w0/training_runs/'
            'design_0000.log')


def model_bref():
    out = {}
    for line in BREF_LOG.read_text().splitlines():
        m = re.match(r'\s*\[b_ref\] (\d{5}) (\d+)\s+(-?[\d.]+) m w\.e\./a', line)
        if m:
            out[(m.group(1), int(m.group(2)))] = float(m.group(3))
    return out


def main():
    acc = pd.read_csv(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else W3 / 'burial_tables'
    F = pd.read_csv(ORI / 'evidence/band_forcing.csv', dtype={'rgi': str})
    band = F.groupby('entity').agg(rgi=('rgi', 'first'), band=('band_model', 'first'))
    n_pts = pd.read_csv(ORI / 'evidence/points.csv').groupby('entity').size()
    A = acc.dropna(subset=['acc_central']).copy()
    unknown = sorted(set(A.entity) - set(band.index))
    if unknown:
        raise SystemExit(f'entities not in the scored set: {unknown}')
    if (A.acc_central <= 0).any():
        raise SystemExit('accumulation must be positive')
    A = A.join(band, on='entity')
    A['n'] = A.entity.map(n_pts)
    A['w'] = A.acc_central * A.n
    g = A.groupby(['rgi', 'band'])
    applied = (g.w.sum() / g.n.sum()).rename('applied')
    A = A.join(applied, on=['rgi', 'band'])
    bref = model_bref()
    A['model_bref'] = [bref.get((r.rgi, int(r.band))) for r in A.itertuples()]
    (out / 'M').mkdir(parents=True, exist_ok=True)
    for rgi, d in applied.groupby(level=0):
        lines = ''.join(f'{int(b):d} {v:.3f}\n' for (_, b), v in d.items())
        text = (f'# measured long-term net accumulation [m w.e./yr], glacier {rgi} (write_burial_tables.py, the '
                f'perfect-burial test of 21 Oct 2026)\n# band elevation [m]  accumulation\n' + lines)
        p = out / 'M' / f'burial_{rgi}.dat'
        if p.exists() and p.read_text() != text:
            raise SystemExit(f'{p} exists with other numbers; use another folder')
        p.write_text(text)
    A.drop(columns=['w']).to_csv(out / 'burial.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(A[['entity', 'rgi', 'band', 'n', 'model_bref', 'acc_central', 'applied', 'confidence']].to_string(index=False))


if __name__ == '__main__':
    main()
