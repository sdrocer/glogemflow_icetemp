"""How well a split of the bands matches the holes' own factors, from the stored pixel factors (no model runs; 22 Oct
2026).

For every scored hole with a week-1 plane (orient_parts5/holes.csv: its own factor F_own, Stage 1 central, dT 0; its
model band; its plane's slope and aspect), the factor of its part under a split (part_factors.SCHEMES; the hole's
part from its plane), computed from the pixels of that band (week3/orient_pixels/pixels, part_factors.py --pixels)
by a chosen statistic: the mean (as in the build), the median, or another quantile. A part without pixels in the
band takes the nearest band of the glacier that has that part, as in the build. Prints the gap to F_own per part of
the hole, also weighted by the hole's scored points, the share of the hole bands' area in each part, and how many
holes' parts come from a neighbouring band.

Usage: python scripts/icetemp_analysis/orientation/pixel_split_test.py [--schemes GNS GNXS30 ...]
           [--stats mean median q60 ...]
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
    if name in ('mean', 'median'):
        return name
    if name.startswith('q'):
        q = float(name[1:]) / 100.0
        return lambda x: x.quantile(q)
    raise ValueError(name)


def lookup(T, rgi, band, part):
    if (rgi, band, part) in T.index:
        return T[(rgi, band, part)], False
    d = T.xs((rgi, part), level=('rgi', 'part'), drop_level=False)
    if d.empty:
        return np.nan, True
    bands = d.index.get_level_values('band').to_numpy()
    return float(d.iloc[int(np.argmin(np.abs(bands - band)))]), True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--schemes', nargs='+', default=['GNS', 'GNESW'])
    ap.add_argument('--stats', nargs='+', default=['mean'])
    a = ap.parse_args()
    H = pd.read_csv(HOLES, dtype={'hole': str, 'rgi': str})
    P = pd.concat([pd.read_csv(PIX / f'pixels_{g}.csv.gz', dtype={'rgi': str},
                               usecols=['rgi', 'band', 'slope', 'aspect', 'F+0']) for g in sorted(H.rgi.unique())],
                  ignore_index=True)
    hole_bands = pd.MultiIndex.from_frame(H[['rgi', 'band']].drop_duplicates())
    w = H.points.to_numpy()
    pd.set_option('display.width', 220)
    rows = []
    for scheme in a.schemes:
        Ps = P.assign(part=pf.classify(np.radians(P.slope.to_numpy()), np.radians(P.aspect.to_numpy()), scheme))
        H['part'] = pf.classify(np.radians(H.slope_deg.to_numpy()), np.radians(H.aspect_deg.to_numpy()), scheme)
        at = Ps.set_index(['rgi', 'band']).loc[lambda d: d.index.isin(hole_bands)]
        share = at.part.value_counts(normalize=True)
        print(f'{scheme}: area share of the hole bands', {p: round(float(share.get(p, 0)), 3) for p in pf.SCHEMES[scheme]},
              '; holes per part', H.part.value_counts().to_dict())
        for stat in a.stats:
            T = Ps.groupby(['rgi', 'band', 'part'])['F+0'].agg(statistic(stat))
            got = [lookup(T, r.rgi, r.band, r.part) for r in H.itertuples()]
            F = pd.Series([g[0] for g in got], index=H.index)
            filled = np.array([g[1] for g in got])
            if stat == 'mean' and scheme in ('GNS', 'GNESW'):
                ref = H.F3 if scheme == 'GNS' else H.F5
                print(f'  mean of the pixels against the build\'s part factors: max difference {(F - ref).abs().max():.4f}')
            e = F - H.F_own
            groups = [('all', np.ones(len(H), bool))] + [(p, (H.part == p).to_numpy()) for p in pf.SCHEMES[scheme]]
            for part, m in groups:
                if not m.any():
                    continue
                rows.append(dict(scheme=scheme, stat=stat, part=part, holes=int(m.sum()), points=int(w[m].sum()),
                                 filled=int(filled[m].sum()), slope=H.slope_deg[m].median(),
                                 F_own=H.F_own[m].median(), F_part=F[m].median(), bias=e[m].mean(),
                                 gap=e[m].abs().mean(), gap_pts=float((e[m].abs() * w[m]).sum() / max(w[m].sum(), 1)),
                                 r=np.corrcoef(F[m], H.F_own[m])[0, 1] if m.sum() > 2 else np.nan))
    R = pd.DataFrame(rows)
    print()
    print(R.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
