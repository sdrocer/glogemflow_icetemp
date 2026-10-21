"""The orientation part (G gentle, N north, S south) of every calibration borehole, for merging the three part runs
(Step 5 plan, step 1c; runner.GloGEMRunner.merge_part_outputs).

A hole takes the part of the EU-DEM 25 m pixel it lies in, with the same slope and aspect as the pixels of
part_factors.py (np.gradient on the EU-DEM, aspect to true north; G below 10 deg, N 270-90 deg, S 90-270 deg), so
holes and pixels are classified alike. For the holes with a week-1 plane (hole_sun.csv, fine DEM within 75 m) the
part from that plane is listed beside it.

Writes OUT/hole_parts.csv (hole, rgi, lat, lon, slope, aspect, part, part_plane).

Usage: PYTHONPATH=src python scripts/icetemp_analysis/orientation/hole_parts.py [OUT]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import exposure_check as ec  # noqa: E402
import forcing_fair_test as ft  # noqa: E402
import part_factors as pf  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler  # noqa: E402


def part_of(slope, aspect):
    return np.where(slope < pf.GENTLE, 'G', np.where((aspect >= np.pi / 2) & (aspect < 1.5 * np.pi), 'S', 'N'))


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else pf.W3 / 'orient_parts'
    out.mkdir(parents=True, exist_ok=True)
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    holes = {}
    for g in dh.calibration_glaciers:
        if g.glacier_id and g.glacier_id not in ft.SEALED:
            for h in g.point_hole:
                holes[str(h)] = g.glacier_id
    B = ts.load_borehole_table().set_index('id')
    eu = ec.EuDem()
    to3035 = Transformer.from_crs('EPSG:4326', 'EPSG:3035', always_xy=True)
    rows = []
    for h, gid in sorted(holes.items()):
        lat, lon = float(B.at[h, 'latitude']), float(B.at[h, 'longitude'])
        x, y = to3035.transform(lon, lat)
        a, x0, y0, cs = eu.block(x - 100, x + 100, y - 100, y + 100)
        dzr, dzc = np.gradient(a, cs)
        r, c = int(round((y0 - y) / cs)), int(round((x - x0) / cs))
        slope = float(np.arctan(np.hypot(dzc[r, c], -dzr[r, c])))
        aspect = float(np.mod(np.arctan2(-dzc[r, c], dzr[r, c]) - ec.north_bearing('EPSG:3035', lon, lat), 2 * np.pi))
        rows.append(dict(hole=h, rgi=gid, lat=lat, lon=lon, slope=np.degrees(slope), aspect=np.degrees(aspect),
                         part=str(part_of(np.array(slope), np.array(aspect)))))
    P = pd.DataFrame(rows)
    H = pd.read_csv(pf.W3 / 'orientation_investigation/evidence/hole_sun.csv', dtype={'hole': str})
    plane = dict(zip(H.hole, part_of(np.radians(H.slope_deg.to_numpy()), np.radians(H.aspect_deg.to_numpy()))))
    P['part_plane'] = P.hole.map(plane)
    P.to_csv(out / 'hole_parts.csv', index=False, float_format='%.4f')
    both = P.dropna(subset=['part_plane'])
    print(f'{len(P)} holes: parts {P.part.value_counts().to_dict()}; with a week-1 plane {len(both)}, '
          f'same part {int((both.part == both.part_plane).sum())}')
    print(pd.crosstab(both.part, both.part_plane, rownames=['EU-DEM 25 m'], colnames=['week-1 plane']).to_string())


if __name__ == '__main__':
    main()
