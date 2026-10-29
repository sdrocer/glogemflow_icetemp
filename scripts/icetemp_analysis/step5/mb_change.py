"""What the free-air mass balance changes for the ice-temperature model (plan of 29 Oct, step 5): per borehole band,
old (2 m-air mass balance, Matthias' calibration) against new (free-air mass balance, its own calibration), from
runs that are already there.

  per glacier  DDFsnow, DDFice, Cprec, T_off, Bw, Ba of both calibration files; the thermal firn line and the ELA,
               mean over the run (firnice_firnline, from the grid runs' gentle part)
  per band     melt     the thermal melt that can refreeze, before the orientation factor [m w.e. per year], mean
                        over the run (firnice_melth of the two melt-logging runs)
               snow_may snow on 1 May [m w.e.], mean over the years (firnice_snowmay, gentle part)
               firn     the share of the year the band is firn (firnice_snowclim firn_frac, mean over months)
               wet      the deepest node the percolating water brought to 0 C [m], mean over the years with a
                        record, and the share of years with one (firnice_wetdepth of the run of the band's part)

The borehole bands are the entities of benchmark_points.csv, each at the nearest GloGEM band of its glacier and the
part of most of its holes (orient_parts5/hole_parts.csv); bands outside the model's glacier take its nearest band
(band_offset, listed in off_model). Writes OUT/mb_change_bands.csv, OUT/mb_change_glaciers.csv
and OUT/mb_change.json (medians of new - old per height range).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/mb_change.py OUT
           --melth OLD_RUN NEW_RUN --grid OLD_BASE NEW_BASE --cal OLD_FILE NEW_FILE
  OLD_BASE / NEW_BASE: grid run names before _oC<P>_w0 (e.g. s5g06_FRHS_r0.85_i1.017_a1.00_fa_fs1.4)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
BAYESCAL = REPO / 'data' / 'bayescal'
FIRNICE = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')
W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
POINTS = REPO / 'data/bayescal/step5_parts5/benchmark/benchmark_points.csv'
HOLE_PARTS = W3 / 'orient_parts5' / 'hole_parts.csv'
RANGES = {'above 3800 m': (3800, np.inf), 'zone 3000-3800 m': (3000, 3800), 'below 3000 m': (-np.inf, 3000)}


def rows(path):
    return [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]


def melth(run, gid):
    """band -> mean yearly thermal melt [m w.e.] over the complete years of the run."""
    r = rows(BAYESCAL / run / FIRNICE / f'firnice_melth_{gid}.dat')
    bands = [float(b) for b in r[0][2:]]
    d = pd.DataFrame([[float(x) for x in l] for l in r[1:]], columns=['year', 'month', *bands])
    full = d.groupby('year').month.count() == 12
    yearly = d[d.year.isin(full.index[full])].groupby('year')[bands].sum()
    return yearly.mean()


def snow_may(run, gid):
    return pd.Series({float(l[0]): np.mean([float(x) for x in l[1:]]) for l in rows(BAYESCAL / run / FIRNICE /
                                                                                     f'firnice_snowmay_{gid}.dat')})


def firn_frac(run, gid):
    d = pd.DataFrame([[float(x) for x in l[:3]] for l in rows(BAYESCAL / run / FIRNICE / f'firnice_snowclim_{gid}.dat')],
                     columns=['elev', 'month', 'firn_frac'])
    return d.groupby('elev').firn_frac.mean()


def wet(run, gid, years):
    """band -> (mean wet depth over the years with a record, share of the run's years with one)."""
    p = BAYESCAL / run / FIRNICE / f'firnice_wetdepth_{gid}.dat'
    d = pd.DataFrame([[float(x) for x in l[:4]] for l in rows(p)], columns=['elev', 'year', 'wet', 'lost'])
    g = d.groupby('elev')
    return g.wet.mean(), g.year.nunique() / years


def firnline(run, gid):
    d = pd.DataFrame([[float(x) for x in l[:4]] for l in rows(BAYESCAL / run / FIRNICE / f'firnice_firnline_{gid}.dat')],
                     columns=['year', 'ela', 'aar', 'firnline_th'])
    return dict(ela=float(d.ela.mean()), firnline_th=float(d.firnline_th.mean()), years=len(d))


def calibration(path):
    d = pd.read_csv(path, sep=r'\s+', dtype={'ID': str})
    return d.set_index(d.ID.str[-5:])[['Ba', 'Bw', 'DDFsnow', 'DDFice', 'Cprec', 'T_off']]


def at(series, band):
    """The value at the nearest band of the series' index."""
    idx = np.asarray(series.index, float)
    return float(series.iloc[int(np.argmin(np.abs(idx - band)))]), float(idx[int(np.argmin(np.abs(idx - band)))])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out')
    ap.add_argument('--melth', nargs=2, required=True, metavar=('OLD_RUN', 'NEW_RUN'))
    ap.add_argument('--grid', nargs=2, required=True, metavar=('OLD_BASE', 'NEW_BASE'))
    ap.add_argument('--cal', nargs=2, required=True, metavar=('OLD_FILE', 'NEW_FILE'))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    P = pd.read_csv(POINTS, dtype={'rgi': str, 'hole': str})
    part_of = dict(pd.read_csv(HOLE_PARTS, dtype={'hole': str})[['hole', 'part']].itertuples(index=False))
    ent = P.groupby('entity').agg(rgi=('rgi', 'first'), band=('band', 'first'), massif=('massif', 'first'),
                                  site=('site', 'first'))
    ent['part'] = [pd.Series([part_of.get(h, 'G') for h in P.hole[P.entity == e]]).mode().iloc[0] for e in ent.index]
    cal = {k: calibration(f) for k, f in zip(('old', 'new'), a.cal)}
    glaciers, bands = [], []
    for gid in sorted(ent.rgi.unique()):
        row = dict(rgi=gid)
        for k, run_m, base in zip(('old', 'new'), a.melth, a.grid):
            for c, v in cal[k].loc[gid].items():
                row[f'{c}_{k}'] = float(v)
            for c, v in firnline(f'{base}_oCG_w0', gid).items():
                row[f'{c}_{k}'] = v
        glaciers.append(row)
        for e in ent.index[ent.rgi == gid]:
            b = float(ent.band[e])
            r = dict(entity=e, rgi=gid, band=b, massif=ent.massif[e], site=ent.site[e], part=ent.part[e])
            for k, run_m, base in zip(('old', 'new'), a.melth, a.grid):
                r[f'melt_{k}'], r['model_band'] = at(melth(run_m, gid), b)
                r[f'snow_may_{k}'], _ = at(snow_may(f'{base}_oCG_w0', gid), b)
                r[f'firn_{k}'], _ = at(firn_frac(f'{base}_oCG_w0', gid), b)
                run_p = f'{base}_oC{ent.part[e]}_w0'
                yrs = firnline(run_p, gid)['years']
                wd, share = wet(run_p, gid, yrs)
                if len(wd) and np.min(np.abs(np.asarray(wd.index, float) - b)) <= 5:
                    r[f'wet_{k}'], _ = at(wd, b)
                    r[f'wet_share_{k}'], _ = at(share, b)
                else:
                    r[f'wet_{k}'], r[f'wet_share_{k}'] = np.nan, 0.0
            bands.append(r)
    B, G = pd.DataFrame(bands), pd.DataFrame(glaciers)
    B['band_offset'] = B.band - B.model_band   # borehole bands outside the model's glacier take its nearest band
    for c in ('melt', 'snow_may', 'firn', 'wet', 'wet_share'):
        B[f'{c}_change'] = B[f'{c}_new'] - B[f'{c}_old']
    for c in ('DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'Ba', 'ela', 'firnline_th'):
        G[f'{c}_change'] = G[f'{c}_new'] - G[f'{c}_old']
    B.to_csv(out / 'mb_change_bands.csv', index=False, float_format='%.4f')
    G.to_csv(out / 'mb_change_glaciers.csv', index=False, float_format='%.4f')
    summ = dict(runs=dict(melth=a.melth, grid=a.grid, cal=a.cal), bands=len(B), glaciers=len(G), ranges={},
                off_model={e: float(o) for e, o in zip(B.entity, B.band_offset) if abs(o) > 5},
                glacier_medians={c: dict(old=float(G[f'{c}_old'].median()), new=float(G[f'{c}_new'].median()),
                                         change=float(G[f'{c}_change'].median()))
                                 for c in ('DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'ela', 'firnline_th')})
    for name, (lo, hi) in RANGES.items():
        d = B[(B.band >= lo) & (B.band < hi)]
        summ['ranges'][name] = dict(bands=len(d), **{c: dict(old=float(d[f'{c}_old'].median()),
                                                             new=float(d[f'{c}_new'].median()),
                                                             change=float(d[f'{c}_change'].median()))
                                                     for c in ('melt', 'snow_may', 'firn', 'wet', 'wet_share')})
    (out / 'mb_change.json').write_text(json.dumps(summ, indent=1))
    pd.set_option('display.width', 200)
    print(G[['rgi'] + [f'{c}_{k}' for c in ('DDFsnow', 'Cprec', 'T_off') for k in ('old', 'new')]].round(2).to_string(index=False))
    for name, v in summ['ranges'].items():
        print(name, v['bands'], {c: f"{v[c]['old']:.2f} -> {v[c]['new']:.2f}" for c in ('melt', 'snow_may', 'firn', 'wet')})
    print(f'-> {out}')


if __name__ == '__main__':
    main()
