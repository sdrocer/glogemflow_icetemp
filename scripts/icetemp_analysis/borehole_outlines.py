"""glenglat boreholes -> the RGI7 outline that contains them, from GloGEM's own glacier masks.

The calibration data matched boreholes to GloGEM glaciers by GLIMS id only, so a hole whose
glenglat GLIMS id is not an RGI7 id was never run (Chessjengletscher lies in 01161, Feegletscher),
and a hole carrying a neighbour's id was modelled on the wrong outline (three Colle Gnifetti holes
labelled 01225 lie in 01230). Rule, per hole, against the rasterised RGI7 outlines GloGEM runs
(rgiv7/glmask/<region>/gl_<id>.grid, UTM 32N):
  inside exactly one mask                          -> that outline            (method 'inside')
  inside several masks                             -> the GLIMS outline if among them, else the first
  outside all, outlines within EDGE_M of the hole  -> the GLIMS outline if among them ('divide'),
                                                      else the nearest one ('nearest')
  nothing within EDGE_M                            -> the GLIMS outline, if any ('glims')
Distances are to the edge of the nearest mask cell. At a divide the masks (40-100 m cells) cannot
say which side a hole is on, so the curator's GLIMS assignment is kept there.

Usage: python scripts/icetemp_analysis/borehole_outlines.py [out_csv]
       default out_csv: data/glenglat_borehole_rgi7.csv (read by icetemp.calibration.data)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

REPO = Path(__file__).resolve().parents[2]
BOREHOLES = REPO / 'glenglat/data/borehole.csv'
LOOKUP = REPO / 'data/rgi7_centraleurope_glims_lookup.csv'
MASKS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/glmask/centraleurope')
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / 'data/glenglat_borehole_rgi7.csv'

EDGE_M = 100.0        # outlines closer than this to a hole outside all masks are candidates
SEARCH_KM = 25.0      # masks of glaciers whose centroid lies within this distance are tested
BOX_MARGIN = 0.2      # deg around the region's glacier centroids: boreholes considered

_masks = {}


def read_mask(gid):
    """(header, 0/1 array, x and y of cell centres) of one glacier mask."""
    if gid not in _masks:
        f = MASKS / f'gl_{gid}.grid'
        if not f.exists():
            _masks[gid] = None
            return None
        hdr = {}
        with open(f) as fh:
            for _ in range(6):
                k, v = fh.readline().split()
                hdr[k.lower()] = float(v)
        a = np.loadtxt(f, skiprows=6) > 0
        rr, cc = np.nonzero(a)
        cs = hdr['cellsize']
        xs = hdr['xllcorner'] + (cc + 0.5) * cs
        ys = hdr['yllcorner'] + (a.shape[0] - 1 - rr + 0.5) * cs
        _masks[gid] = (hdr, a, xs, ys)
    return _masks[gid]


def locate(gid, x, y):
    """(inside, distance to the nearest cell edge [m]) of point x, y for glacier gid's mask."""
    m = read_mask(gid)
    if m is None or len(m[2]) == 0:
        return False, np.inf
    hdr, a, xs, ys = m
    cs = hdr['cellsize']
    c = int(np.floor((x - hdr['xllcorner']) / cs))
    r = int(a.shape[0] - 1 - np.floor((y - hdr['yllcorner']) / cs))
    if 0 <= r < a.shape[0] and 0 <= c < a.shape[1] and a[r, c]:
        return True, 0.0
    dd = np.hypot(np.maximum(np.abs(xs - x) - cs / 2, 0), np.maximum(np.abs(ys - y) - cs / 2, 0))
    return False, float(dd.min())


def assign(row, look, tr):
    x, y = tr.transform(row.longitude, row.latitude)
    near = np.hypot((look.cenlat - row.latitude) * 111.0,
                    (look.cenlon - row.longitude) * 111.0 * np.cos(np.deg2rad(row.latitude)))
    cand = look.gid[near < SEARCH_KM].tolist()
    by_glims = look.gid[look.glims_id == row.glims_id]
    g = by_glims.iloc[0] if len(by_glims) else None
    res = {gid: locate(gid, x, y) for gid in cand}
    inside = [gid for gid, (ins, _) in res.items() if ins]
    if len(inside) == 1:
        return inside[0], 'inside', 0.0, g, ''
    if len(inside) > 1:
        pick = g if g in inside else inside[0]
        return pick, 'inside', 0.0, g, ','.join(inside)
    close = sorted((d, gid) for gid, (_, d) in res.items() if d <= EDGE_M)
    if close:
        dist = dict((gid, d) for d, gid in close)
        if g in dist:
            return g, 'divide', dist[g], g, ','.join(gid for _, gid in close)
        return close[0][1], 'nearest', close[0][0], g, ','.join(gid for _, gid in close)
    return g, 'glims' if g else 'none', np.nan, g, ''


def main():
    look = pd.read_csv(LOOKUP)
    look['gid'] = look.rgi7_glogem.str.split('.').str[-1]
    bh = pd.read_csv(BOREHOLES)
    box = ((bh.latitude >= look.cenlat.min() - BOX_MARGIN) & (bh.latitude <= look.cenlat.max() + BOX_MARGIN)
           & (bh.longitude >= look.cenlon.min() - BOX_MARGIN) & (bh.longitude <= look.cenlon.max() + BOX_MARGIN))
    bh = bh[box & bh.latitude.notna() & bh.longitude.notna()]
    tr = Transformer.from_crs('EPSG:4326', 'EPSG:32632', always_xy=True)
    rows = []
    for row in bh.itertuples():
        gid, method, dist, g, others = assign(row, look, tr)
        rows.append(dict(borehole_id=row.id, glacier_name=row.glacier_name, glims_id=row.glims_id,
                         rgi7_by_glims=g, rgi7=gid, method=method,
                         edge_distance_m=None if dist is None or not np.isfinite(dist) else round(dist, 1),
                         candidates=others))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    moved = out[out.rgi7.notna() & (out.rgi7 != out.rgi7_by_glims)]
    print(f'{len(out)} boreholes in the region box -> {OUT}')
    print(out.method.value_counts().to_string())
    print(f'outline differs from the GLIMS lookup for {len(moved)} holes:')
    print(moved[['borehole_id', 'glacier_name', 'rgi7_by_glims', 'rgi7', 'method', 'edge_distance_m']].to_string(index=False))


if __name__ == '__main__':
    main()
