"""How well a part factor matches the holes' own factors, from the stored pixel factors (no model runs; 22 Oct 2026).

For every scored hole with a week-1 plane (orient_parts5/holes.csv: its own factor F_own, Stage 1 central, dT 0; its
model band; its parts under three and five parts), the factor of its part computed from the pixels of that band
(week3/orient_pixels/pixels, part_factors.py --pixels) by a chosen statistic: the mean (as in the build), the median,
or another quantile. A part without pixels in the band takes the nearest band of the glacier that has that part, as
in the build. Prints the gap to F_own per part of the hole, also weighted by the hole's scored points.

Usage: python scripts/icetemp_analysis/orientation/pixel_split_test.py [--stats mean median q60 ...]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import part_factors as pf  # noqa: E402

PIX = pf.W3 / 'orient_pixels/pixels'
HOLES = pf.W3 / 'orient_parts5/holes.csv'


def statistic(name):
    if name == 'mean':
        return 'mean'
    if name == 'median':
        return 'median'
    if name.startswith('q'):
        q = float(name[1:]) / 100.0
        return lambda x: x.quantile(q)
    raise ValueError(name)


def part_table(P, scheme, stat):
    P = P.assign(part=pf.classify(np.radians(P.slope.to_numpy()), np.radians(P.aspect.to_numpy()), scheme))
    return P.groupby(['rgi', 'band', 'part'])['F+0'].agg(statistic(stat))


def lookup(T, rgi, band, part):
    if (rgi, band, part) in T.index:
        return T[(rgi, band, part)], band
    d = T.xs((rgi, part), level=('rgi', 'part'), drop_level=False)
    if d.empty:
        return np.nan, np.nan
    bands = d.index.get_level_values('band').to_numpy()
    k = int(np.argmin(np.abs(bands - band)))
    return float(d.iloc[k]), bands[k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stats', nargs='+', default=['mean', 'median'])
    a = ap.parse_args()
    H = pd.read_csv(HOLES, dtype={'hole': str, 'rgi': str})
    P = pd.concat([pd.read_csv(PIX / f'pixels_{g}.csv.gz', dtype={'rgi': str},
                               usecols=['rgi', 'band', 'slope', 'aspect', 'F+0']) for g in sorted(H.rgi.unique())],
                  ignore_index=True)
    pd.set_option('display.width', 200)
    for scheme, col, ref in (('GNS', 'part3', 'F3'), ('GNESW', 'part5', 'F5')):
        out = []
        for stat in a.stats:
            T = part_table(P, scheme, stat)
            got = [lookup(T, r.rgi, r.band, getattr(r, col)) for r in H.itertuples()]
            H[f'F_{stat}'] = [g[0] for g in got]
            if stat == 'mean':
                dev = (H.F_mean - H[ref]).abs().max()
                print(f'{scheme}: mean of the pixels against the build\'s part factors: max difference {dev:.4f}')
            e = H[f'F_{stat}'] - H.F_own
            w = H.points.to_numpy()
            for part, d in [('all', H.index)] + [(p, H.index[H[col] == p]) for p in pf.SCHEMES[scheme]]:
                ee, ww = e.loc[d], w[H.index.get_indexer(d)]
                out.append(dict(scheme=scheme, stat=stat, part=part, holes=len(d), points=int(ww.sum()),
                                F_own=H.F_own.loc[d].median(), F_part=H[f'F_{stat}'].loc[d].median(),
                                bias=ee.mean(), gap=ee.abs().mean(),
                                gap_pts=float((ee.abs() * ww).sum() / max(ww.sum(), 1)),
                                r=np.corrcoef(H[f'F_{stat}'].loc[d], H.F_own.loc[d])[0, 1] if len(d) > 2 else np.nan))
        R = pd.DataFrame(out)
        print(R.pivot_table(index='part', columns='stat', values=['F_part', 'bias', 'gap', 'gap_pts', 'r'],
                            sort=False).round(3).to_string())
        print(R[R.stat == a.stats[0]].set_index('part')[['holes', 'points', 'F_own']].round(3).T.to_string(), '\n')


if __name__ == '__main__':
    main()
