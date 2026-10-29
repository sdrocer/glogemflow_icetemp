"""The yardstick of 29 Oct: model versions scored on the same 1,135 readings below 15 m, each massif predicted from
the other six (plan 2026-10-29_freeair_mb_icetemp_plan.md). Reported per version:

  M1, zone_M1, above_M1   the benchmark's error (outlines equal, >= 3 points), overall, 3000-3800 m, >= 3800 m
  rmse, r2, bias          per reading, overall and per height range (r2 = 1 - SS_res / SS_tot)
  level, shape            per profile (one hole in one year): its mean, and its deviations from that mean
                          (profiles with >= 3 readings over >= 10 m); grad = change with depth, C per 100 m
  classes                 bands in the right regime class per height range (one tolerance per band, from its
                          first hole: approximate, chapter 16 counts 73 of 105 for step5 where this gives 72)
  bars                    bar 1 = o5, bar 2 = o4 (M1, above_M1, zone_M1); 'as good' = at most 0.1 C above the bar
                          overall and above 3800 m, zone at most 0.2 C above

Versions are LABEL=FILE:COLUMN; every file must list the readings in the order of the reference points file, whose
hole and year columns define the profiles. --check asserts the M1 of the four versions known on 29 Oct.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/yardstick.py [--add LABEL=FILE:COL]
           [--out FILE] [--check]
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step5_common as sc  # noqa: E402
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402

WEEK3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
REF = sc.REPO / 'data/bayescal/step5_parts5/benchmark/benchmark_points.csv'
VERSIONS = {
    'o4 (2 m air, three parts)': f'{WEEK3}/freeair/score/points.csv:BASE',
    'o5 (2 m air, five parts)': f'{WEEK3}/freeair/score_fb1/points.csv:BASE',
    'fb1 (free air, five parts)': f'{WEEK3}/freeair/score_fb1/points.csv:FA',
    'option A grid choice': f'{REF}:grid',
    'option A Step 5': f'{REF}:step5',
    'v1.0 in-sample (0.85 / 1.017)': f'{REF}:G@0.85',
    'line in elevation': f'{REF}:elev_line',
}
KNOWN_M1 = {'o4 (2 m air, three parts)': 1.855, 'o5 (2 m air, five parts)': 2.126,
            'fb1 (free air, five parts)': 2.791, 'option A grid choice': 2.831}
BAR1, BAR2 = 'o5 (2 m air, five parts)', 'o4 (2 m air, three parts)'
RANGES = {'above 3800 m': (3800, np.inf), 'zone 3000-3800 m': (3000, 3800), 'below 3000 m': (-np.inf, 3000)}


def m1(d, c):
    g = ((d[c] - d.obs) ** 2).groupby(d.rgi)
    mse = g.mean()[g.size() >= 3]
    return float(np.sqrt(mse.mean()))


def errors(d, c):
    e = d[c] - d.obs
    return dict(n=len(d), rmse=float(np.sqrt((e ** 2).mean())), bias=float(e.mean()),
                r2=float(1 - (e ** 2).sum() / ((d.obs - d.obs.mean()) ** 2).sum()))


def profiles(d, c):
    g = d.groupby('prof')
    lev = g.agg(obs=('obs', 'mean'), mod=(c, 'mean'), n=('obs', 'size'), span=('depth', lambda x: x.max() - x.min()))
    keep = lev.index[(lev.n >= 3) & (lev.span >= 10)]
    sh = d[d.prof.isin(keep)]
    od, md = sh.obs - sh.groupby('prof').obs.transform('mean'), sh[c] - sh.groupby('prof')[c].transform('mean')
    grad = np.array([(np.polyfit(q.depth, q.obs, 1)[0], np.polyfit(q.depth, q[c], 1)[0]) for _, q in sh.groupby('prof')]) * 100
    return dict(profiles=len(lev), level_rmse=float(np.sqrt(((lev['mod'] - lev.obs) ** 2).mean())),
                level_r=float(np.corrcoef(lev.obs, lev['mod'])[0, 1]), shape_profiles=int(len(keep)),
                shape_rmse=float(np.sqrt(((md - od) ** 2).mean())), shape_r=float(np.corrcoef(od, md)[0, 1]),
                grad_measured=float(np.median(grad[:, 0])), grad_model=float(np.median(grad[:, 1])))


def classes(d, c, tol, rng_of):
    cls = {e: ft.classify_entity(q[c].to_numpy(), q.depth.to_numpy(), 'hole', tol[e]) for e, q in d.groupby('entity')}
    obs = {e: ft.classify_entity(q.obs.to_numpy(), q.depth.to_numpy(), 'hole', tol[e]) for e, q in d.groupby('entity')}
    out = {r: [0, 0] for r in RANGES}
    for e in cls:
        out[rng_of[e]][1] += 1
        out[rng_of[e]][0] += int(cls[e] == obs[e])
    return {r: f'{a} of {b}' for r, (a, b) in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--add', action='append', default=[], help='LABEL=FILE:COLUMN, scored beside the known versions')
    ap.add_argument('--out', default=str(sc.REPO / 'data/bayescal/step5_parts5/benchmark/yardstick.json'))
    ap.add_argument('--check', action='store_true', help='assert the M1 of the versions known on 29 Oct')
    a = ap.parse_args()
    versions = dict(VERSIONS)
    for spec in a.add:
        label, src = spec.split('=', 1)
        versions[label] = src
    S = pd.read_csv(REF, dtype={'rgi': str, 'hole': str})
    S['prof'] = S.entity + '|' + S.hole + '|' + S.year.astype(str)
    cache = {}
    for label, src in versions.items():
        path, col = src.rsplit(':', 1)
        if path not in cache:
            cache[path] = pd.read_csv(path, dtype={'rgi': str})
            F = cache[path]
            if len(F) != len(S) or not (F.entity.to_numpy() == S.entity.to_numpy()).all() \
                    or not np.allclose(F.depth.to_numpy(), S.depth.to_numpy()) \
                    or not np.allclose(F.obs.to_numpy(), S.obs.to_numpy()):
                raise SystemExit(f'{path}: readings differ from {REF}')
        S[label] = cache[path][col].to_numpy()
    if S[list(versions)].isna().any().any():
        raise SystemExit('a version has missing values')
    bh = ts.load_borehole_table()
    tol = {e: ts.borehole_tolerance(q.hole.iloc[0], bh) for e, q in S.groupby('entity')}
    band = S.groupby('entity').band.first()
    rng_of = {e: next(r for r, (lo, hi) in RANGES.items() if lo <= b < hi) for e, b in band.items()}
    res = dict(commit=subprocess.run(['git', '-C', str(sc.REPO), 'rev-parse', 'HEAD'], capture_output=True,
                                     text=True).stdout.strip(),
               points=len(S), sources=versions, versions={})
    for label in versions:
        v = dict(M1=m1(S, label), zone_M1=m1(S[S.band.between(*ft.ZONE)], label),
                 above_M1=m1(S[S.band >= 3800], label), all=errors(S, label), classes=classes(S, label, tol, rng_of))
        for r, (lo, hi) in RANGES.items():
            d = S[(S.band >= lo) & (S.band < hi)]
            v[r] = dict(**errors(d, label), **profiles(d, label))
        res['versions'][label] = v
    bars = {}
    for name, ref in [('bar 1', BAR1), ('bar 2', BAR2)]:
        b = res['versions'][ref]
        bars[name] = {label: dict(M1=v['M1'] - b['M1'], above_M1=v['above_M1'] - b['above_M1'],
                                  zone_M1=v['zone_M1'] - b['zone_M1'],
                                  as_good=bool(v['M1'] - b['M1'] <= 0.1 and v['above_M1'] - b['above_M1'] <= 0.1
                                               and v['zone_M1'] - b['zone_M1'] <= 0.2))
                      for label, v in res['versions'].items()}
    res['bars'] = bars
    if a.check:
        for label, want in KNOWN_M1.items():
            if abs(res['versions'][label]['M1'] - want) > 5e-4:
                raise SystemExit(f'{label}: M1 {res["versions"][label]["M1"]:.4f}, expected {want}')
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(f'{"version":32s} {"M1":>5s} {"above":>5s} {"zone":>5s} {"rmse":>5s} {"r2":>5s}  bar 1  bar 2')
    for label, v in res['versions'].items():
        print(f'{label:32s} {v["M1"]:5.2f} {v["above_M1"]:5.2f} {v["zone_M1"]:5.2f} {v["all"]["rmse"]:5.2f} '
              f'{v["all"]["r2"]:5.2f}  {"yes" if bars["bar 1"][label]["as_good"] else "no ":5s}  '
              f'{"yes" if bars["bar 2"][label]["as_good"] else "no"}')
    print(f'-> {a.out}')


if __name__ == '__main__':
    main()
