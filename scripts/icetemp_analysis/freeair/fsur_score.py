"""The f_sur test (27 Oct 2026; book ch. 14): can f_ref and f_sur together put the high sites right under the free air?

Runs: free air, five orientation parts (factors rebuilt for the free air) merged per hole, f_ins 1.017, f_adv 1.0, on
the grid f_sur 0 / 1.5 / 3.0 / 4.5 K x f_ref 0.40 / 0.55 / 0.70 / 0.85 / 1.00 (f_sur 0 = the base fb1; the rest tag
fs1). The pair (f_sur, f_ref) is chosen leaving one massif out (lowest error on the other six); the base is the same
choice over f_ref alone at f_sur 0.

The reading, fixed before the runs (revised after the max-effort review, before any fs1 run: per band and site
instead of pooled means; f_ref 0.40 added; 30.0 m readings counted as deep):
- A band at or above 3800 m fits if its mean error is within +-1.7 C (the floor between neighbouring boreholes) in each
  depth class where it has at least 3 readings: shallow (15-20 m) and deep (30-200 m). A site fits if more than half
  of its high bands fit (sites from the site table).
- f_ref + f_sur fix the high sites if (1) more high sites fit than with the base, and at least half of them fit;
  (2) f_ref is below 1.0 in at least 4 of the 7 folds and f_sur below 4.5 in at least 4 of the 7 folds (a choice at a
  grid edge means the range is too small, not a fit).
Reported, not part of the reading: the f_sur each high site needs at each f_ref (interpolated, deep and shallow;
none where the site responds by less than 0.1 C per K of f_sur), the
Mont Blanc minus Monte Rosa gap of the deep error, the pooled errors above 3800 m, the error on left-out massifs
(overall, above 3800 m, zone), the regime classes at every grid point and out of sample.
Checked before scoring: every fs1 part run on GloGEM 61ada09 without local changes, and f_sur applied (the last
month's surface T_sfc equals fb1's minus f_sur on every band where fb1's is below 0).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/freeair/fsur_score.py [OUT_DIR]
           [--dry]   (--dry: f_sur 0 only, checks the scorer against fb1)
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
SUB = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')
FSUR = (0.0, 1.5, 3.0, 4.5)
FREF = (0.40, 0.55, 0.70, 0.85, 1.00)
PARTS = 'GNESW'
GLOGEM = '61ada09'
HIGH, FLOOR, MIN_READ = 3800, 1.7, 3
MIN_SLOPE = 0.1   # C per K of f_sur: below it a site does not respond and gets no needed f_sur


def run_name(fs, fr, part='P'):
    if fs == 0:
        return f'fb1_FRHS_r{fr:.2f}_i1.017_a1.00_fa_oC{part}_w0'
    return f'fs1_FRHS_r{fr:.2f}_i1.017_a1.00_fa_fs{fs:.1f}_oC{part}_w0'


def read_tsfc(run):
    out = {}
    for f in (ft.BAYESCAL / run / SUB).glob('firnice_temp_validation_*.dat'):
        lines = f.read_text().splitlines()
        cols = next(ln for ln in lines if ln.startswith('# elev')).lstrip('#').split()
        v = np.array([ln.split() for ln in lines if ln.strip() and not ln.startswith('#')
                      and len(ln.split()) == len(cols)], float)
        for e, t in zip(v[:, cols.index('elev')], v[:, cols.index('T_sfc')]):
            out[(f.stem.split('_')[-1], e)] = t
    return out


def check_runs(fsur):
    n_checked = 0
    for fs in [x for x in fsur if x > 0]:
        for fr in FREF:
            for p in PARTS:
                run = run_name(fs, fr, p)
                prov = (ft.BAYESCAL / run / 'provenance.txt').read_text().split()
                if not prov[1].startswith(GLOGEM) or 'none' not in prov:
                    raise SystemExit(f'{run}: GloGEM {prov[1][:7]} or local changes (want {GLOGEM}, none)')
                a, b = read_tsfc(run_name(0, fr, p)), read_tsfc(run)
                keys = [k for k in a if a[k] < 0 and k in b]
                bad = [k for k in keys if abs(b[k] - (a[k] - fs)) > 0.002]
                if bad or not keys:
                    raise SystemExit(f'{run}: f_sur not applied on {len(bad)} of {len(keys)} bands, e.g. {bad[:3]}')
                n_checked += len(keys)
    return n_checked


def band_errors(Q, col):
    """Per band at or above 3800 m: mean error in the shallow (<= 20 m) and deep (>= 30 m) class, NaN with fewer than
    MIN_READ readings."""
    rows = []
    for e, d in Q[Q.band >= HIGH].groupby('entity'):
        err = d[col] - d.obs
        sh, dp = err[d.depth <= 20], err[d.depth >= 30]
        rows.append(dict(entity=e, site=d.site.iloc[0], massif=d.massif.iloc[0],
                         shallow=sh.mean() if len(sh) >= MIN_READ else np.nan,
                         deep=dp.mean() if len(dp) >= MIN_READ else np.nan))
    B = pd.DataFrame(rows)
    B = B[B[['shallow', 'deep']].notna().any(axis=1)].copy()
    B['fits'] = ((B.shallow.abs() <= FLOOR) | B.shallow.isna()) & ((B.deep.abs() <= FLOOR) | B.deep.isna())
    return B


def sites_fitting(B):
    s = B.groupby('site').fits.agg(['sum', 'size'])
    return sorted(s.index[s['sum'] > s['size'] / 2]), len(s)


def massif_gap(B):
    m = B.groupby('massif').deep.mean()
    return float(m.get('Mont Blanc', np.nan) - m.get('Monte Rosa', np.nan))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'freeair' / 'fsur'))
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    fsur = (0.0,) if a.dry else FSUR
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    grid = [(fs, fr) for fs in fsur for fr in FREF]
    runs = {ft.col(f'S{fs:.1f}', fr): run_name(fs, fr) for fs, fr in grid}
    missing = [r for r in runs.values() if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} outputs missing, e.g. {sorted(missing)[:3]}')
    n_checked = check_runs(fsur)

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    zone, hi = P.band.between(*ft.ZONE), P.band >= HIGH
    obs_cls = ft.entity_classes(P, 'obs', E)
    names = sorted(P.entity.unique())

    def summary(cols):
        pred, chosen = ft.lomo(P, cols)
        Q = P.assign(oos=pred)
        B = band_errors(Q, 'oos')
        fit, n_sites = sites_fitting(B)
        cls = ft.entity_classes(Q, 'oos', E)
        e = Q.oos - Q.obs
        return dict(chosen=chosen, M1=ft.m1(Q, 'oos'), above_M1=ft.m1(Q[hi], 'oos'), zone_M1=ft.m1(Q[zone], 'oos'),
                    pooled_15_20=float(e[hi & (P.depth <= 20)].mean()), pooled_30_200=float(e[hi & (P.depth >= 30)].mean()),
                    high_sites=n_sites, sites_fit=fit, n_sites_fit=len(fit), mont_blanc_minus_monte_rosa=massif_gap(B),
                    counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, cls, names).items()}), B

    base, Bb = summary({fr: ft.col('S0.0', fr) for fr in FREF})
    test, Bt = summary({(fs, fr): ft.col(f'S{fs:.1f}', fr) for fs, fr in grid})
    ch = list(test['chosen'].values())
    test['folds_f_ref_below_1'] = sum(fr < 1.0 for fs, fr in ch)
    test['folds_f_sur_below_top'] = sum(fs < max(fsur) for fs, fr in ch)
    test['edge_choices'] = sorted({f'{m}: f_sur {fs}, f_ref {fr}' for m, (fs, fr) in test['chosen'].items()
                                   if fr in (FREF[0], FREF[-1]) or fs == max(fsur)})
    test['fixes_high_sites'] = bool(test['n_sites_fit'] > base['n_sites_fit'] and
                                    test['n_sites_fit'] >= test['high_sites'] / 2 and
                                    test['folds_f_ref_below_1'] >= 4 and test['folds_f_sur_below_top'] >= 4)
    Bb.to_csv(out / 'bands_high_base.csv', index=False, float_format='%.3f')
    Bt.to_csv(out / 'bands_high_test.csv', index=False, float_format='%.3f')

    pts, need = [], []
    for fs, fr in grid:
        c = ft.col(f'S{fs:.1f}', fr)
        B = band_errors(P.assign(m=P[c]), 'm')
        cls = ft.entity_classes(P, c, E)
        pts.append(dict(f_sur=fs, f_ref=fr, M1_in_sample=ft.m1(P, c), n_sites_fit=len(sites_fitting(B)[0]),
                        gap=massif_gap(B), **{f'cls_{k}': f'{x}/{n}' for k, (x, n) in
                                              ft.class_counts(obs_cls, cls, names).items()}))
    if len(fsur) > 1:
        for fr in FREF:
            site_err = {fs: band_errors(P.assign(m=P[ft.col(f'S{fs:.1f}', fr)]), 'm').groupby('site')[['shallow', 'deep']].mean()
                        for fs in fsur}
            for site in site_err[0.0].index:
                row = dict(site=site, f_ref=fr)
                for k in ('shallow', 'deep'):
                    y = np.array([site_err[fs].loc[site, k] for fs in fsur])
                    if np.isfinite(y).all():
                        slope, icpt = np.polyfit(fsur, y, 1)
                        # a site that barely responds to f_sur (temperate, or melt-dominated) has no needed value
                        row[f'f_sur_needed_{k}'] = float(-icpt / slope) if abs(slope) >= MIN_SLOPE else np.nan
                        row[f'slope_{k}'] = float(slope)
                need.append(row)
    pd.DataFrame(pts).to_csv(out / 'grid_points.csv', index=False, float_format='%.3f')
    pd.DataFrame(need).to_csv(out / 'f_sur_needed.csv', index=False, float_format='%.2f')
    res = dict(points=len(P), grid=grid, runs=runs, runs_checked_bands=n_checked,
               base={k: v for k, v in base.items() if k != 'chosen'} | dict(chosen=base['chosen']),
               test={k: v for k, v in test.items() if k != 'chosen'} |
                    dict(chosen={m: dict(f_sur=fs, f_ref=fr) for m, (fs, fr) in test['chosen'].items()}))
    (out / 'fsur_score.json').write_text(json.dumps(res, indent=1, default=str))
    for lab, r in (('base (f_sur 0)', base), ('f_ref + f_sur', test)):
        print(f"{lab}: high sites fitting {r['n_sites_fit']} of {r['high_sites']} {r['sites_fit']}; "
              f"left-out massifs {r['M1']:.2f} (above 3800 m {r['above_M1']:.2f}); pooled above 3800 m "
              f"{r['pooled_15_20']:+.2f} / {r['pooled_30_200']:+.2f}; Mont Blanc - Monte Rosa deep {r['mont_blanc_minus_monte_rosa']:+.2f}")
    print('chosen:', test['chosen'], '| edges:', test['edge_choices'], '| fixes the high sites:', test['fixes_high_sites'])


if __name__ == '__main__':
    main()
