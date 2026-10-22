"""Three parts or five? A check without model runs (21-22 Oct 2026, Janosch: "what if we also include west and east").

For every scored borehole with a week-1 plane: its own orientation factor (Stage 1, central, dT 0: the hole's plane
and horizon; orient_tables/factors.csv), and the factor of the part it falls in, at its model band, under the
three-part split (gentle / north 270-90 / south 90-270; orient_parts_v2) and under the five-part split (gentle /
north, east, south, west as 90-deg sectors; orient_parts5). The better split gives each hole a factor closer to its
own. Also, per band, how far the west part's factor lies from the east part's, and from north and south.

Writes OUT/holes.csv and prints the summary.
Usage: python scripts/icetemp_analysis/orientation/parts_scheme_compare.py [OUT]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import part_factors as pf  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
ORI = W3 / 'orientation_investigation'


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'orient_parts5'
    H = pd.read_csv(ORI / 'evidence/hole_sun.csv', dtype={'hole': str, 'rgi': str})
    own = pd.read_csv(W3 / 'orient_tables/factors.csv', dtype={'hole': str})
    own = own[(own.strength == 'C') & (own.dT == 0.0)].set_index('hole').F
    bm = pd.read_csv(ORI / 'evidence/band_forcing.csv', dtype={'rgi': str}).groupby('entity').band_model.first()
    pts = pd.read_csv(ORI / 'evidence/points.csv', dtype={'hole': str}).groupby('hole').size()
    P3 = pd.read_csv(W3 / 'orient_parts_v2/parts.csv', dtype={'rgi': str}).set_index(['rgi', 'band', 'part'])['F+0']
    P5r = pd.read_csv(W3 / 'orient_parts5/parts.csv', dtype={'rgi': str})
    P5 = P5r.set_index(['rgi', 'band', 'part'])['F+0']
    sl, asp = np.radians(H.slope_deg.to_numpy()), np.radians(H.aspect_deg.to_numpy())
    H['part3'], H['part5'] = pf.classify(sl, asp, 'GNS'), pf.classify(sl, asp, 'GNESW')
    H['band'] = H.entity.map(bm).astype(float)
    H['F_own'] = H.hole.map(own)
    H['F3'] = [P3.get((r.rgi, r.band, r.part3), np.nan) for r in H.itertuples()]
    H['F5'] = [P5.get((r.rgi, r.band, r.part5), np.nan) for r in H.itertuples()]
    H['points'] = H.hole.map(pts).fillna(0)
    H = H.dropna(subset=['F_own', 'F3', 'F5'])
    H[['hole', 'rgi', 'entity', 'band', 'slope_deg', 'aspect_deg', 'part3', 'part5', 'F_own', 'F3', 'F5',
       'points']].to_csv(out / 'holes.csv', index=False, float_format='%.4f')

    def stats(d):
        w = d.points.to_numpy()
        e3, e5 = (d.F3 - d.F_own).abs(), (d.F5 - d.F_own).abs()
        return pd.Series(dict(holes=len(d), points=int(w.sum()), r3=np.corrcoef(d.F3, d.F_own)[0, 1] if len(d) > 2 else np.nan,
                              r5=np.corrcoef(d.F5, d.F_own)[0, 1] if len(d) > 2 else np.nan,
                              gap3=float(e3.mean()), gap5=float(e5.mean()),
                              gap3_pts=float((e3 * w).sum() / max(w.sum(), 1)),
                              gap5_pts=float((e5 * w).sum() / max(w.sum(), 1))))
    pd.set_option('display.width', 200)
    print('all scored holes with a plane:')
    print(stats(H).round(3).to_string())
    print('\nby five-part sector of the hole:')
    print(H.groupby('part5').apply(stats).round(3).to_string())
    # per band: west against east, north and south (bands where the parts have their own pixels)
    R = P5r[P5r.n_pix >= 5].pivot_table(index=['rgi', 'band'], columns='part', values='F+0')
    R = R[R.index.get_level_values('band') >= 3000]
    for a, b in (('W', 'E'), ('W', 'N'), ('W', 'S'), ('E', 'N')):
        d = (R[a] - R[b]).dropna()
        print(f'bands >= 3000 m with both {a} and {b}: {len(d)}; {a} minus {b}: median {d.median():+.3f}, '
              f'5-95 % {d.quantile(.05):+.3f} to {d.quantile(.95):+.3f}')
    share = P5r[P5r.band >= 3000].groupby('part').n_pix.sum()
    print('pixel shares above 3000 m:', (share / share.sum()).round(3).to_dict())


if __name__ == '__main__':
    main()
