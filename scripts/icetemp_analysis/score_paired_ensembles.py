"""Paired scores of two ensembles run on one design: what changes from ensemble a to ensemble b at
the same parameter sets and the same borehole points (Step 4: alps_r7 against alps_s2).

Usage: score_paired_ensembles.py --a <label>=<source> --b <label>=<source> [--out <file.md>]
                                 [--design <design.csv>]

A source is an ensemble tag prefix (reads data/bayescal/<tag>_w<k>/training_runs/
design_XXXX.output.pkl for every worker k) or a directory of re-matched design_XXXX.pkl files
(e.g. tmp_icetemp_trees/rescored/alps_r7_step0). Each file maps entity -> (depths, model values),
the entity's observation points first. The design is design.csv of the worker folders; a directory
without its own design.csv needs --design, and every design.csv found must agree with it. To score
one ensemble alone, pass it twice under two labels. Outputs written before icetemp 2f603aa
(8 Oct 2026) were matched at the end of the month: pass their re-matched directory (as
alps_r7_step0 for alps_r7), never the tag; the band check reports such outputs as not reproduced.

Compared are the design points present in both ensembles, at the common points: (entity, point)
finite in every one of them in both ensembles. Residual = model - observed. The report (Markdown,
printed and written to --out): the model band each entity was compared with (tag sources: the raw
profiles archived for one design point, re-matched file by file); pooled RMSE and bias; the change
per design point (b - a); borehole groups, firn label and depth, for the best run and the median
over runs; parameter strengths (per point, a linear fit over the design points with f_ref >= 0.1,
in prior sd units) and how alike their patterns are; the best 5 runs; every design point.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/score_paired_ensembles.py \
    --a r7=/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/rescored/alps_r7_step0 \
    --b s2=alps_s2 --design data/bayescal/alps_r7_w0/design.csv
"""
import argparse
import pickle
import re
import sys
import tarfile
import tempfile
import warnings
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from icetemp.calibration.data import DataHandler, _elevation_band
from icetemp.calibration.priors import Priors
from icetemp.calibration.runner import PROFILE_ARCHIVE_DIR, match_model_to_obs
from icetemp.io import read_profile_file

REPO = Path(__file__).resolve().parents[2]
BAYESCAL = REPO / 'data' / 'bayescal'
PARAMS = ['refreeze_frac', 'insul_scale', 'advection_scale']
SHORT = ['f_ref', 'f_ins', 'f_adv']
PAIRS = [(0, 1), (0, 2), (1, 2)]
# the 7 test glaciers (ensemble_sweep.ENSEMBLE)
TEST_GLACIERS = ['Grenzgletscher', 'Glacier des Bossons', 'Glacier de Tête Rousse', 'Hintereisferner',
                 'Gornergletscher', 'Hohsaasgletscher', 'Vadret dal Corvatsch']
MONTE_ROSA = ['01225', '01230']
# groups split on the entity's 10 m band centre (...5 m), so no entity sits on a split
HIGH = 3000            # m; tongues below, firn above
CG_FROM = 4400         # m; Colle Gnifetti saddle (glenglat CG); below: Grenz basin GG, Lysjoch LJ, Seserjoch SJ
CG_DEEP = 40.0         # m; Colle Gnifetti deeper than this: climate before 1942, accepted
F_REF_MIN = 0.1        # design points below are left out of the parameter fit
MIN_FIT = 8            # design points needed for the 4-term fit
OUT_RES = 0.001        # degC, resolution of GloGEM's profile files (f8.3)
NOISE_X = 4            # a cosine needs both slope sets above this many times their rounding noise
SHALLOW = 15.0         # m; strengths of firn and ice points also split here (seasonal wave, winter snow)
DEPTH_BINS = [(-np.inf, 5, '0-5 m'), (5, 15, '5-15 m'), (15, 30, '15-30 m'), (30, np.inf, '> 30 m')]
N_BEST = 5


def load_source(spec, flag):
    """label, runs {design index: output dict}, origin {design index: worker folder}, designs
    (design.csv files), where, bad (unreadable files), is_dir."""
    label, sep, src = spec.partition('=')
    if not (sep and label and src):
        raise SystemExit(f'{flag} wants <label>=<source>, got {spec!r}')
    path = Path(src)
    if path.is_dir():
        files = [(f, None) for f in sorted(path.glob('design_*.pkl'))]
        designs = [path / 'design.csv'] if (path / 'design.csv').is_file() else []
        where = f'directory {path.resolve()}'
    else:
        workers = sorted((d for d in BAYESCAL.glob(f'{src}_w*')
                          if re.fullmatch(re.escape(src) + r'_w\d+', d.name)),
                         key=lambda d: int(d.name.rsplit('_w', 1)[1]))
        if not workers:
            raise SystemExit(f'{src}: neither a directory nor an ensemble {BAYESCAL}/{src}_w<k>')
        files = [(f, d.name) for d in workers for f in sorted((d / 'training_runs').glob('design_*.output.pkl'))]
        designs = [d / 'design.csv' for d in workers if (d / 'design.csv').is_file()]
        where = f'ensemble {BAYESCAL}/{src}_w<k>, k = {", ".join(d.name.rsplit("_w", 1)[1] for d in workers)}'
    runs, origin, seen, bad = {}, {}, {}, []
    for f, worker in files:
        m = re.fullmatch(r'design_(\d+)(?:\.output)?\.pkl', f.name)
        if m is None:
            continue
        i = int(m.group(1))
        if i in seen:
            raise SystemExit(f'{label}: design point {i} found twice ({seen[i]}, {f})')
        seen[i] = f
        try:
            with open(f, 'rb') as fh:
                runs[i] = pickle.load(fh)
            origin[i] = worker
        except Exception as err:          # e.g. a file a running worker is still writing
            bad.append(f'{f} ({type(err).__name__})')
    return SimpleNamespace(label=label, runs=runs, origin=origin, designs=designs, where=where, bad=bad,
                           is_dir=path.is_dir())


def read_design(path):
    d = pd.read_csv(path)
    missing = [p for p in PARAMS if p not in d.columns]
    if missing:
        raise SystemExit(f'{path}: no column {missing}')
    return d[PARAMS].to_numpy(float)


def model_matrix(runs, idx, ents, label):
    """Model value at every observation point (rows) for every design point (columns); NaN where
    the entity or the depth is missing."""
    cols = []
    for i in idx:
        col = []
        for g in ents:
            v = np.full(g.n_obs, np.nan)
            o = runs[i].get(g.glacier_name)
            if o is not None:
                d, T = np.asarray(o[0], float), np.asarray(o[1], float)
                same = len(T) == len(d) >= g.n_obs and np.allclose(d[:g.n_obs], g.depths, rtol=0, atol=1e-6)
                if not same:
                    raise SystemExit(f'{label} design {i}: {g.glacier_name} was matched at other depths '
                                     f'than the observation points loaded now')
                v = T[:g.n_obs]
            col.append(v)
        cols.append(np.concatenate(col))
    return np.column_stack(cols)


def identical_runs(M, idx):
    """Design point pairs with identical model values everywhere: one output was parsed from
    another run's profiles."""
    return [(idx[j], idx[k]) for j in range(M.shape[1]) for k in range(j + 1, M.shape[1])
            if np.array_equal(M[:, j], M[:, k], equal_nan=True)]


def same_output(x, y):
    return all(np.shape(u) == np.shape(v) and np.array_equal(np.asarray(u, float), np.asarray(v, float),
                                                              equal_nan=True) for u, v in zip(x, y))


def band_check(src, idx, ents):
    """Re-match the raw profiles archived for the first paired design point that has an archive
    (runner.archive_profiles) with every equally near band file: which band each entity was
    compared with, and whether its stored values are reproduced. None without an archive."""
    cands = [(i, PROFILE_ARCHIVE_DIR / src.origin[i] / f'design_{i:04d}.tar.gz') for i in idx if src.origin.get(i)]
    hit = next(((i, p) for i, p in cands if p.is_file()), None)
    if hit is None:
        return None
    i, tgz = hit
    prof = {}
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(tgz) as tar:
            keep = [m for m in tar.getmembers() if re.fullmatch(r'temp_ID\d+_\w+\.dat', m.name)]
            tar.extractall(tmp, members=keep, filter='data')
        for f in sorted(Path(tmp).glob('temp_ID*.dat')):
            elev, _, df = read_profile_file(f)
            prof.setdefault(f.stem.rsplit('_', 1)[1], []).append((elev, df))
    used, lost = {}, []
    for g in ents:
        stored, cand = src.runs[i].get(g.glacier_name), prof.get(g.glacier_id, [])
        if stored is None or not cand:
            continue
        dmin = min(abs(e - g.elevation) for e, _ in cand)
        near = [(e, df) for e, df in cand if abs(e - g.elevation) == dmin]
        hits = sorted({e for e, df in near if not df.empty and same_output(match_model_to_obs(df, g), stored)})
        if hits:
            used[g.glacier_name] = (hits, sorted({e for e, _ in near}), dmin)
        else:
            lost.append(g.glacier_name)
    return SimpleNamespace(design=i, archive=tgz, used=used, lost=lost)


def rmse(R):
    return np.sqrt(np.mean(R ** 2, axis=0))


def cosine(u, v):
    with np.errstate(invalid='ignore', divide='ignore'):
        return float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))


def num(x, nd=3, sign=False):
    if x is None or not np.isfinite(x):
        return '-'
    return f'{x:+.{nd}f}' if sign else f'{x:.{nd}f}'


def md_table(headers, rows):
    """Markdown table padded for the terminal; first column left, the others right."""
    cells = [[str(c) for c in headers]] + [[str(c) for c in r] for r in rows]
    w = [max(len(r[k]) for r in cells) for k in range(len(headers))]
    fmt = lambda r: '| ' + ' | '.join(c.ljust(w[k]) if k == 0 else c.rjust(w[k]) for k, c in enumerate(r)) + ' |'
    sep = '|' + '|'.join(':' + '-' * (n + 1) if k == 0 else '-' * (n + 1) + ':' for k, n in enumerate(w)) + '|'
    return [fmt(cells[0]), sep] + [fmt(r) for r in cells[1:]]


MR_FIRN = f'Monte Rosa high firn ({HIGH}-{CG_FROM} m)'
CG = f'Colle Gnifetti (>= {CG_FROM} m)'


def borehole_groups(P):
    r, z, gl = P.rgi, P.band, P.glacier
    mr = r.isin(MONTE_ROSA)
    cg = mr & (z >= CG_FROM)
    main = [
        (f'Hintereisferner tongue, ice (< {HIGH} m)', r.eq('03116') & (z < HIGH)),
        (f'Hintereisferner firn (>= {HIGH} m)', r.eq('03116') & (z >= HIGH)),
        (f'Grenz tongue (< {HIGH} m)', mr & (z < HIGH) & gl.eq('Grenzgletscher')),
        ('Gorner tongue', gl.eq('Gornergletscher')),
        ('Tête Rousse', r.eq('00777')),
        ('Hohsaas', r.eq('01112')),
        (MR_FIRN, mr & (z >= HIGH) & (z < CG_FROM) & ~gl.eq('Gornergletscher')),
        (CG, cg),
        ('Mont Blanc / Bossons', r.eq('00773')),
        ('Corvatsch', r.eq('02103')),
    ]
    other = ~np.logical_or.reduce([m.to_numpy() for _, m in main])
    out = [('all common points', np.ones(len(P), bool))] + [(n, m.to_numpy()) for n, m in main]
    k = [n for n, _ in out].index(CG) + 1
    out[k:k] = [(f'  of it 0-{CG_DEEP:.0f} m', (cg & (P.depth <= CG_DEEP)).to_numpy()),
                (f'  of it deeper than {CG_DEEP:.0f} m (accepted)', (cg & (P.depth > CG_DEEP)).to_numpy())]
    if other.any():
        out.append(('other', other))
    return out


def site_counts(boreholes, P, m):
    """glenglat label prefixes of the boreholes in the 10 m bands of the entities under m."""
    keys = set(zip(P.glacier[m], P.band[m]))
    sel = [(n, b) in keys for n, b in zip(boreholes.glacier_name, boreholes.band)]
    pre = boreholes.label[sel].astype(str).str.extract(r'^([A-Z]{2})\d', expand=False).fillna('other')
    return ', '.join(f'{k} {v}' for k, v in pre.value_counts().sort_index().items())


def label_and_depth(P):
    out = [(f'{lab} points', P.label.eq(lab).to_numpy()) for lab in ('firn', 'ice', 'unknown')]
    out += [(name, ((P.depth > lo) & (P.depth <= hi)).to_numpy()) for lo, hi, name in DEPTH_BINS]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--a', required=True, metavar='LABEL=SOURCE', help='reference ensemble')
    ap.add_argument('--b', required=True, metavar='LABEL=SOURCE', help='changed ensemble')
    ap.add_argument('--out', help='Markdown file for the report')
    ap.add_argument('--design', help='design.csv; needed for a directory source without one')
    args = ap.parse_args()

    sa, sb = load_source(args.a, '--a'), load_source(args.b, '--b')
    la, lb = sa.label, sb.label
    if la == lb:
        raise SystemExit('the two labels must differ')
    for s in (sa, sb):
        if s.is_dir and not s.designs and not args.design:
            raise SystemExit(f'{s.label}: {s.where} has no design.csv: pass --design (the design its runs used)')
    design_files = ([Path(args.design)] if args.design else []) + sa.designs + sb.designs
    if not design_files:
        raise SystemExit('no design.csv in either source: pass --design')
    design = read_design(design_files[0])
    for f in design_files[1:]:
        d = read_design(f)
        if d.shape != design.shape or not np.allclose(d, design, rtol=0, atol=1e-9):
            raise SystemExit(f'{f} differs from {design_files[0]}: not one shared design')
    outside = {s.label: sorted(i for i in s.runs if i >= len(design)) for s in (sa, sb)}
    idx = sorted(set(sa.runs) & set(sb.runs))
    if not idx:
        raise SystemExit(f'no design point in both ensembles ({la}: {len(sa.runs)}, {lb}: {len(sb.runs)})')
    if idx[-1] >= len(design):
        raise SystemExit(f'design points {[i for i in idx if i >= len(design)]} are not in the '
                         f'{len(design)}-point design')
    for lab, out in outside.items():
        if out:
            print(f'WARNING: {lab} has design points outside the {len(design)}-point design (ignored): '
                  f'{", ".join(map(str, out))}', file=sys.stderr)
    n_des = len(idx)
    theta = design[idx]

    print('loading the borehole data (about 2 min) ...', file=sys.stderr, flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        dh = DataHandler(region='CentralEurope', include_estimated=True, include_unassessed=True)
        dh.load()
    known = {g.glacier_name: g for g in dh.calibration_glaciers}
    names = {n for i in idx for o in (sa.runs[i], sb.runs[i]) for n in o}
    ents = sorted((known[n] for n in names if n in known and known[n].base_glacier_name in TEST_GLACIERS),
                  key=lambda g: (g.base_glacier_name, g.elevation))
    not_in_data = sorted(n for n in names if n not in known)
    not_test = sorted(n for n in names if n in known and known[n].base_glacier_name not in TEST_GLACIERS)
    not_run = sorted(g.glacier_name for g in dh.calibration_glaciers
                     if g.base_glacier_name in TEST_GLACIERS and g.glacier_name not in names)
    if not ents:
        raise SystemExit('no entity of the 7 test glaciers in the outputs')

    pts = pd.DataFrame([dict(entity=g.glacier_name, glacier=g.base_glacier_name, rgi=g.glacier_id,
                             elev=g.elevation, band=g.elevation_band, depth=float(z), label=str(lab),
                             obs=float(t))
                        for g in ents for z, t, lab in zip(g.depths, g.T_obs, g.firn_label)])
    A = model_matrix(sa.runs, idx, ents, la)
    B = model_matrix(sb.runs, idx, ents, lb)
    fin_a, fin_b = np.isfinite(A), np.isfinite(B)
    common = fin_a.all(axis=1) & fin_b.all(axis=1)
    if not common.any():
        raise SystemExit('no point is finite in every paired design point of both ensembles')
    P = pts[common].reset_index(drop=True)
    y = P.obs.to_numpy()[:, None]
    RA, RB = A[common] - y, B[common] - y
    rm_a, rm_b = rmse(RA), rmse(RB)
    bi_a, bi_b = RA.mean(axis=0), RB.mean(axis=0)
    ja, jb = int(np.argmin(rm_a)), int(np.argmin(rm_b))
    wa, wb = int(np.argmax(rm_a)), int(np.argmax(rm_b))
    dups = {la: identical_runs(A, idx), lb: identical_runs(B, idx)}
    checks = {s.label: (None if s.is_dir else band_check(s, idx, ents)) for s in (sa, sb)}
    n_common = P.entity.value_counts()
    n_drop = pts.entity[~common].value_counts()

    L = []
    L += [f'# Paired ensemble scores: {la} (a) against {lb} (b)', '',
          f'- a = {la}: {sa.where}; {len(sa.runs)} design points read',
          f'- b = {lb}: {sb.where}; {len(sb.runs)} design points read',
          f'- design: {design_files[0]} ({len(design)} points)'
          + (f'; the other {len(design_files) - 1} design files agree' if len(design_files) > 1 else ''),
          f'- paired design points (in both): {n_des} of {len(design)}']
    for s in (sa, sb):
        miss = sorted(set(range(len(design))) - set(s.runs))
        if miss:
            L.append(f'  - not in {s.label}: {", ".join(map(str, miss))}')
        if outside[s.label]:
            L.append(f'  - **in {s.label} but not in the design (ignored): '
                     f'{", ".join(map(str, outside[s.label]))}**')
        if s.bad:
            L.append(f'  - unreadable in {s.label} (skipped): {"; ".join(s.bad)}')
    if n_des < len(design):
        L += ['', f'**PARTIAL: only {n_des} of {len(design)} design points are in both ensembles. '
                  'Every number below will change when the rest are in.**']
    for lab, d in dups.items():
        if d:
            L += ['', f'**WARNING: in {lab}, design points with identical model values at every point (one '
                      f'output parsed from another run\'s profiles?): {", ".join(f"{i}={j}" for i, j in d)}**']
    for lab, c in checks.items():
        if c is not None and c.lost:
            some = ', '.join(c.lost[:6]) + (f' and {len(c.lost) - 6} more' if len(c.lost) > 6 else '')
            L += ['', f'**WARNING: {lab} design point {c.design}: the stored values of {len(c.lost)} entities '
                      'are not reproduced from its archived profiles by today\'s matcher (outputs from an older '
                      f'comparison, or another run\'s profiles): {some}**']
    L += ['', f'- borehole data: DataHandler(CentralEurope, estimated and unassessed included), '
              f'{len(ents)} entities of the 7 test glaciers in the outputs',
          f'- common points (finite in all {n_des} paired design points of both ensembles): '
          f'{int(common.sum())} of {len(pts)}; dropped {int((~common).sum())}']
    if not_run:
        L.append(f'- entities of the 7 test glaciers not in the outputs: {", ".join(not_run)}')
    if not_test:
        L.append(f'- in the outputs but not one of the 7 test glaciers (left out): {", ".join(not_test)}')
    if not_in_data:
        L.append(f'- in the outputs but not in the borehole data (left out): {", ".join(not_in_data)}')
    L += ['- residual = model - observed (°C); RMSE and bias over the common points',
          f'- written {datetime.now():%Y-%m-%d %H:%M}', '']

    # model band per entity
    L += ['## Model band per entity', '',
          'The matcher (runner.parse_training_output) compares an entity with the model column of the band '
          'nearest its elevation (that of its first borehole). On a band edge two bands are equally near, and '
          'only one is used. Checked by re-matching the archived raw profiles of one design point against '
          'every equally near band:', '']
    off_band = {}
    for s in (sa, sb):
        c = checks[s.label]
        if s.is_dir:
            L.append(f'- {s.label}: not checked (a directory source has no profile archive)')
            continue
        if c is None:
            L.append(f'- {s.label}: not checked (no profile archive under {PROFILE_ARCHIVE_DIR})')
            continue
        L.append(f'- {s.label}, design point {c.design} ({c.archive}): stored values reproduced for '
                 f'{len(c.used)} of {len(c.used) + len(c.lost)} entities')
        off_band[s.label] = {e: u for e, u in c.used.items() if u[0] != [known[e].elevation_band]}
    rows = []
    for lab, ob in off_band.items():
        for e, (hits, near, dmin) in sorted(ob.items()):
            g = known[e]
            why = (f'on a band edge: {" and ".join(map(str, near))} equally near' if g.elevation_band in near
                   else f'no band of its own; nearest {dmin:.0f} m away')
            rows.append([e, lab, g.glacier_id, f'{g.elevation:.1f}', g.elevation_band, '/'.join(map(str, hits)),
                         why, int(n_common.get(e, 0)), int(n_drop.get(e, 0))])
    if rows:
        L += ['', 'Entities compared with a band other than their own 10 m band:', '']
        L += md_table(['entity', 'checked', 'RGI7', 'elevation (m)', 'own band', 'band used', 'why',
                       'common points', 'dropped points'], rows)
    elif off_band:
        L += ['', 'Every checked entity was compared with its own 10 m band.']
    L += ['']

    if (~common).any():
        L += ['## Points dropped', '']
        rows = []
        for r in np.flatnonzero(~common):
            e = pts.entity[r]
            row = [e, pts.rgi[r], f'{pts.depth[r]:.1f}', pts.label[r], f'{int((~fin_a[r]).sum())} of {n_des}',
                   f'{int((~fin_b[r]).sum())} of {n_des}']
            for lab, ob in off_band.items():
                row.append(f'{"/".join(map(str, ob[e][0]))} (own {known[e].elevation_band})' if e in ob else
                           'own' if e in checks[lab].used else '?')
            rows.append(row)
        L += md_table(['entity', 'RGI7', 'depth (m)', 'label', f'not finite in {la}', f'not finite in {lb}']
                      + [f'band used ({lab})' for lab in off_band], rows)
        L += ['', 'Not finite: the depth lies outside the model column used (above its first node or below its '
                  'deepest resolved node), or the entity is missing from that run. An entity read from a '
                  'neighbouring band can lose depths that its own band\'s column reaches (each band\'s column '
                  'runs as deep as the deepest reading in that band).', '']

    # 1. pooled
    L += ['## 1. Pooled RMSE and bias', '', 'Best = lowest RMSE among the paired design points.', '']
    rows = [['best RMSE (design point)', f'{num(rm_a[ja])} ({idx[ja]})', f'{num(rm_b[jb])} ({idx[jb]})',
             num(rm_b[jb] - rm_a[ja], sign=True)],
            ['bias at the best run', num(bi_a[ja], sign=True), num(bi_b[jb], sign=True), ''],
            ['median RMSE', num(np.median(rm_a)), num(np.median(rm_b)),
             num(np.median(rm_b) - np.median(rm_a), sign=True)],
            ['median bias', num(np.median(bi_a), sign=True), num(np.median(bi_b), sign=True),
             num(np.median(bi_b) - np.median(bi_a), sign=True)],
            ['worst RMSE (design point)', f'{num(rm_a[wa])} ({idx[wa]})', f'{num(rm_b[wb])} ({idx[wb]})',
             num(rm_b[wb] - rm_a[wa], sign=True)],
            ['bias at the worst run', num(bi_a[wa], sign=True), num(bi_b[wb], sign=True), ''],
            ['bias range', f'{num(bi_a.min(), sign=True)} .. {num(bi_a.max(), sign=True)}',
             f'{num(bi_b.min(), sign=True)} .. {num(bi_b.max(), sign=True)}', '']]
    L += md_table(['°C', la, lb, f'{lb} - {la}'], rows) + ['']

    # 2. paired change
    dr, db = rm_b - rm_a, bi_b - bi_a
    dm = B[common] - A[common]
    r, j = np.unravel_index(int(np.argmax(np.abs(dm))), dm.shape)
    L += ['## 2. Change per design point (b - a)', '',
          f'Each design point scored in both ensembles; {lb} - {la}, over the {n_des} paired design points.', '']
    rows = [['RMSE', num(np.median(dr), sign=True), f'{num(dr.min(), sign=True)} .. {num(dr.max(), sign=True)}',
             f'{int((dr < 0).sum())} of {n_des} lower in {lb}'],
            ['bias', num(np.median(db), sign=True), f'{num(db.min(), sign=True)} .. {num(db.max(), sign=True)}',
             f'{int((np.abs(bi_b) < np.abs(bi_a)).sum())} of {n_des} closer to 0 in {lb}']]
    L += md_table(['change (°C)', 'median', 'range', 'count'], rows)
    if dm[r, j] == 0:
        L += ['', 'Largest change of one model value at a common point: 0 (identical).', '']
    else:
        L += ['', f'Largest change of one model value at a common point: {dm[r, j]:+.4g} °C '
                  f'({P.entity[r]}, {P.depth[r]:.1f} m, design point {idx[j]}).', '']

    def breakdown(masks, with_rgi):
        best, med = [], []
        for name, m in masks:
            if not m.any():
                continue
            ga, gb = rmse(RA[m]), rmse(RB[m])
            ba, bb = RA[m].mean(axis=0), RB[m].mean(axis=0)
            d = gb - ga
            head = [name] + (['/'.join(sorted({str(x) for x in P.rgi[m]}))] if with_rgi else []) + [int(m.sum())]
            best.append(head + [num(ga[ja], 2), num(ba[ja], 2, True), num(gb[jb], 2), num(bb[jb], 2, True)])
            med.append(head + [num(np.median(ga), 2), num(np.median(gb), 2),
                               f'{num(np.median(d), 2, True)} ({num(d.min(), 2, True)} .. '
                               f'{num(d.max(), 2, True)})',
                               num(np.median(ba), 2, True), num(np.median(bb), 2, True),
                               num(np.median(bb - ba), 2, True), num(np.abs(dm[m]).max(), 4)])
        h = [''] + (['RGI7'] if with_rgi else []) + ['points']
        return (md_table(h + [f'{la} RMSE', f'{la} bias', f'{lb} RMSE', f'{lb} bias'], best),
                md_table(h + [f'{la} RMSE', f'{lb} RMSE', f'RMSE change: median (range)', f'{la} bias',
                              f'{lb} bias', 'bias change', 'largest change of one value'], med))

    groups = borehole_groups(P)
    for title, masks, with_rgi in (('3. Borehole groups', groups, True),
                                   ('4. Firn label and depth', label_and_depth(P), False)):
        t_best, t_med = breakdown(masks, with_rgi)
        L += [f'## {title}', '']
        if with_rgi:
            L += [f'Grouped by the entity\'s 10 m band (the @ in its name): tongues below {HIGH} m; on '
                  f'Monte Rosa the saddle of Colle Gnifetti from {CG_FROM} m, the firn of the Grenz basin, '
                  f'Lysjoch and Seserjoch below it. Colle Gnifetti deeper than {CG_DEEP:.0f} m reflects the '
                  'climate before 1942, which the run does not have (accepted, shown apart).', '']
        L += [f'Best run of each ensemble ({la}: design {idx[ja]}, {lb}: design {idx[jb]}):', '']
        L += t_best + ['', f'Median over the {n_des} paired design points; the change is {lb} - {la} '
                           'per design point, then its median and range; the last column is the largest '
                           f'|{lb} - {la}| of one model value over the points and design points:', '']
        L += t_med + ['']
        if with_rgi:
            bh = pd.read_csv(Path(dh.glenglat_dir) / 'borehole.csv',
                             usecols=['glacier_name', 'elevation', 'label'])
            bh = bh[bh.elevation.notna()]
            bh = bh.assign(band=[_elevation_band(float(e)) for e in bh.elevation])
            gm = dict(groups)
            L += ['glenglat boreholes in the bands of the Monte Rosa firn groups, by label (CG Colle Gnifetti, '
                  'GG Grenzgletscher, LJ Lysjoch, SJ Seserjoch): '
                  + '; '.join(f'{n}: {site_counts(bh, P, gm[n]) or "none"}' for n in (MR_FIRN, CG)) + '.', '']

    # 5. parameter strengths
    L += ['## 5. Parameter strengths', '']
    keep = theta[:, 0] >= F_REF_MIN
    if keep.sum() < MIN_FIT:
        L += [f'Skipped: {int(keep.sum())} paired design points with f_ref >= {F_REF_MIN}, '
              f'{MIN_FIT} needed for the fit.', '']
    else:
        pri = Priors()
        sd = np.array([pri.refreeze_frac.std(), pri.insul_scale.std(), pri.advection_scale.std()])
        X = np.column_stack([np.ones(int(keep.sum())), theta[keep] / sd])
        S = {lab: np.linalg.lstsq(X, M[common][:, keep].T, rcond=None)[0][1:].T for lab, M in ((la, A), (lb, B))}
        # slope noise from rounding alone: each value off by up to OUT_RES/2 (sd OUT_RES/sqrt 12)
        floor = NOISE_X * OUT_RES / np.sqrt(12) * np.sqrt(np.diag(np.linalg.inv(X.T @ X))[1:])
        moved = {lab: np.ptp(M[common][:, keep], axis=1) > OUT_RES + 1e-9 for lab, M in ((la, A), (lb, B))}
        left = [str(i) for i, k in zip(idx, keep) if not k]
        L += [f'Per point, the model value over the {int(keep.sum())} paired design points with f_ref >= '
              f'{F_REF_MIN} (left out: {", ".join(left) or "none"}) fitted by a plane in the three parameters, '
              f'each in prior standard deviations (f_ref {sd[0]:.4f}, f_ins {sd[1]:.4f}, f_adv {sd[2]:.4f}). '
              'Strength = the slope, °C per prior sd: median of |slope| over the points. Unmoved = points whose '
              f'model value changes by at most {OUT_RES} °C (the output resolution) over these design points, '
              'e.g. ice held at the melting point: their slopes are 0.', '']

        def strength_rows(masks):
            rows = []
            for name, m in masks:
                if m.any():
                    rows.append([name, int(m.sum())] + [
                        x for lab in (la, lb) for x in [int((m & ~moved[lab]).sum())]
                        + [num(np.median(np.abs(S[lab][m, k]))) for k in range(3)]])
            return md_table(['', 'points'] + [f'{lab} {p}' for lab in (la, lb) for p in ['unmoved'] + SHORT],
                            rows)

        shallow = (P.depth <= SHALLOW).to_numpy()
        masks = []
        for lab in ('firn', 'ice'):
            m = P.label.eq(lab).to_numpy()
            masks += [(f'{lab} points', m), (f'  of it 0-{SHALLOW:.0f} m', m & shallow),
                      (f'  of it deeper than {SHALLOW:.0f} m', m & ~shallow)]
        masks += [('unknown points', P.label.eq('unknown').to_numpy()), ('all points', np.ones(len(P), bool))]
        L += ['By firn label and depth:', ''] + strength_rows(masks)
        L += ['', 'By borehole group (as in section 3):', ''] + strength_rows(groups[1:]) + ['']

        rows_of = {e: r for e, r in P.groupby('entity', sort=False).indices.items() if len(r) >= 3}
        cls = {}
        for e, r_ in rows_of.items():
            share = P.label.iloc[r_].value_counts(normalize=True)
            cls[e] = share.index[0] if share.iloc[0] > 0.5 else 'mixed'
        rms = lambda u: np.sqrt(np.mean(u ** 2))
        rows = []
        for p, q in PAIRS:
            for kind in ('all', 'firn', 'ice', 'unknown', 'mixed'):
                sel = [r_ for e, r_ in rows_of.items() if kind == 'all' or cls[e] == kind]
                if not sel:
                    continue
                row = [f'{SHORT[p]} ~ {SHORT[q]}' if kind == 'all' else '', f'{kind} entities', len(sel)]
                for lab in (la, lb):
                    c = np.array([cosine(S[lab][r_, p], S[lab][r_, q]) for r_ in sel if moved[lab][r_].sum() >= 3
                                  and rms(S[lab][r_, p]) >= floor[p] and rms(S[lab][r_, q]) >= floor[q]])
                    row += [len(c)] + ([num(np.median(c), 2, True), num(np.median(np.abs(c)), 2)] if len(c) else
                                       ['-', '-'])
                rows.append(row)
        L += ['How alike the patterns are, per entity: the cosine between the slopes of two parameters over the '
              'entity\'s points, i.e. whether they bend that one profile into the same shape (+1 / -1: the same '
              'shape, that profile alone cannot tell them apart; 0: unrelated). n = entities with >= 3 common '
              'points; one is used for a pair where >= 3 of its points move and both parameters move them by at '
              f'least {NOISE_X} times the slope that rounding the output to {OUT_RES} °C gives by itself (RMS '
              f'over its points; here f_ref {floor[0]:.1e}, f_ins {floor[1]:.1e}, f_adv {floor[2]:.1e} °C per '
              'prior sd); below that, its shape is noise. Entities are classed by the label of most of their '
              'points. Median over the entities used:', '']
        L += md_table(['', 'entities', 'n', f'{la} n used', f'{la} cos', f'{la} abs cos', f'{lb} n used',
                       f'{lb} cos', f'{lb} abs cos'], rows)
        rows = []
        for name, m in (('all points', np.ones(len(P), bool)), ('firn points', P.label.eq('firn').to_numpy()),
                        ('ice points', P.label.eq('ice').to_numpy())):
            if m.any():
                rows.append([name, int(m.sum())] + [num(cosine(S[lab][m, p], S[lab][m, q]), 2, True)
                                                     for lab in (la, lb) for p, q in PAIRS])
        L += ['', 'The same cosine over pooled points, across sites at once: it also counts where each '
                  'parameter acts (f_ins on ice, where f_ref does nothing), which a single profile cannot show. '
                  'Both cosines compare shapes, not sizes: a parameter that barely acts (f_ref on ice, see the '
                  'strengths) can share the other\'s shape and still not matter. Neither is the test of '
                  'whether a calibration can tell the parameters apart; that is Step 5a.', '']
        L += md_table(['', 'points'] + [f'{lab} {SHORT[p]} ~ {SHORT[q]}' for lab in (la, lb) for p, q in PAIRS],
                      rows) + ['']

    # 6. best runs
    L += [f'## 6. The best {N_BEST} runs', '',
          'Design ranges: ' + ', '.join(f'{s} {design[:, k].min():.2f}-{design[:, k].max():.2f}'
                                        for k, s in enumerate(SHORT)) + '.', '']
    for lab, rm, bi, other, olab in ((la, rm_a, bi_a, rm_b, lb), (lb, rm_b, bi_b, rm_a, la)):
        rows = [[n + 1, idx[j]] + [num(theta[j, k], 2) for k in range(3)] +
                [num(rm[j]), num(bi[j], sign=True), num(other[j])]
                for n, j in enumerate(np.argsort(rm, kind='stable')[:N_BEST])]
        L += [f'{lab}:', '']
        L += md_table(['rank', 'design point'] + SHORT + ['RMSE', 'bias', f'RMSE in {olab}'], rows) + ['']

    # appendix
    L += ['## Every paired design point', '']
    rows = [[i] + [num(theta[j, k], 2) for k in range(3)] +
            [num(rm_a[j]), num(rm_b[j]), num(dr[j], sign=True), num(bi_a[j], sign=True), num(bi_b[j], sign=True),
             num(db[j], sign=True)] for j, i in enumerate(idx)]
    L += md_table(['design point'] + SHORT + [f'{la} RMSE', f'{lb} RMSE', 'change', f'{la} bias', f'{lb} bias',
                                              'change'], rows)

    text = '\n'.join(L) + '\n'
    print(text)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding='utf-8')
        print(f'written to {out}', file=sys.stderr)


if __name__ == '__main__':
    main()
