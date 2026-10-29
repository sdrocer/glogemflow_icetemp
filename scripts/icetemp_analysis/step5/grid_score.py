"""The grid runs of a Step 5 config scored alone, before any fit (plan of 29 Oct, step 4): the f_ref grid 0.03-1.00 at
f_ins 1.017 and the starting values, each massif predicted from the other six, with the benchmark's own functions
(benchmark.py: frame, add_run, the part checks of grid_runs / start_run, forcing_fair_test.lomo, run_measures,
class_measures, frozen_beds).

Writes OUT/grid_points.csv (the benchmark's points, columns G@<f_ref>, start and grid; in benchmark_points.csv's
order, so yardstick.py can add it) and OUT/grid_score.json:
  grid_chosen   f_ref chosen for each left-out massif
  best_own      each massif's own best f_ref on the grid (in-sample, rule 1a)
  methods       grid and start: M1 overall, zone, above 3800 m; classes; frozen beds

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/grid_score.py CONFIG [--out DIR]
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark as bm  # noqa: E402
import step5_common as sc  # noqa: E402
from icetemp.calibration.config import CalibrationConfig  # noqa: E402
from icetemp.calibration.pipeline import CalibrationPipeline  # noqa: E402

ft, ts = bm.ft, bm.ts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config')
    ap.add_argument('--out', default=None, help='default: <output_dir>/grid')
    a = ap.parse_args()
    head = subprocess.run(['git', '-C', str(sc.REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    pl = CalibrationPipeline(CalibrationConfig.from_yaml(a.config))
    pl.load_data()
    out = Path(a.out) if a.out else pl.config.output_path / 'grid'
    out.mkdir(parents=True, exist_ok=True)
    glaciers = pl._calibration_glaciers()
    E = ft.scored_entities(glaciers, ts.load_borehole_table(), pl.config.fit_exclusions)
    site_of = pl.site_of(glaciers)
    P = bm.frame(E, glaciers, site_of)
    gr = bm.grid_runs(pl.config)
    cols = {fr: ft.col('G', fr) for fr in bm.GRID}
    for fr, run in gr.items():
        bm.add_run(P, cols[fr], run)
    start = bm.start_run(pl.config)
    bm.add_run(P, 'start', start)
    P['grid'], chosen = ft.lomo(P, cols)
    methods = ['grid', 'start']
    bad = {c: sorted(P.entity[~np.isfinite(P[c].to_numpy(float))].unique())
           for c in methods + list(cols.values()) if not np.isfinite(P[c].to_numpy(float)).all()}
    if bad:
        raise SystemExit(f'predictions missing: {bad}')
    names = sorted(P.entity.unique())
    obs_cls = bm.classes_of(P, 'obs', E)
    best_own = {m: float(min(bm.GRID, key=lambda fr: ft.m1(P[P.massif == m], cols[fr])))
                for m in bm.massifs_with_points(P)}
    res = dict(commit=head, config=str(Path(a.config).resolve()), grid_runs=gr, start_run=start, points=len(P),
               entities=len(names), grid_chosen=chosen, best_own=best_own, methods={})
    for c in methods:
        res['methods'][c] = dict(**bm.run_measures(P, c),
                                 classes=bm.class_measures(obs_cls, bm.classes_of(P, c, E), names, site_of))
    hp = dict(pd.read_csv(bm.HOLE_PARTS, dtype={'hole': str})[['hole', 'part']].itertuples(index=False))
    res['methods']['grid']['frozen_beds'], _ = bm.frozen_beds(E, chosen, {fr: gr[fr][:-len('_oCP_w0')] for fr in bm.GRID},
                                                              hp)
    res['methods']['start']['frozen_beds'], _ = bm.frozen_beds(E, {m: 'start' for m in bm.massifs_with_points(P)},
                                                               {'start': start[:-len('_oCP_w0')]}, hp)
    P.to_csv(out / 'grid_points.csv', index=False, float_format='%.3f')
    (out / 'grid_score.json').write_text(json.dumps(res, indent=1, default=float))
    for c in methods:
        m = res['methods'][c]
        print(f'{c:6s} M1 {m["M1"]:.3f}  zone {m["zone_M1"]:.3f}  above {m["above_M1"]:.3f}  '
              f'bands right {m["classes"]["bands_right"]}/{m["classes"]["bands"]}  frozen beds {m["frozen_beds"]}')
    print('grid chosen', chosen)
    print('best own  ', best_own)
    print(f'-> {out}')


if __name__ == '__main__':
    main()
