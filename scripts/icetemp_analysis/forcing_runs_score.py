"""Score the thermal forcing runs of the plan of 2026-10-12 against today's run, with the checks fixed
before the runs.

Runs (data/bayescal/<tag>_w0, prior-mean parameters, all 24 calibration glaciers): today =
ce25_pred_w0 (before the switches), R (new code, defaults), F (no t_offset in the thermal forcing),
H (Herron & Langway 1980 + latent heat per m water), FH (both), P (H + water stops at the model's
own 830 kg/m3 depth; diagnostic).

Checks: R byte-identical to today; the mass-balance files of F, H, FH unchanged; F leaves the
outlines with t_offset 0.00 byte-identical; Colle Gnifetti 01230 minus 01225 at 10-40 m below
1 C (F, FH); the correlation of deep bias with t_offset across outlines no longer significant
(F, FH; permutation test); the deep RMSE of the glaciers outside Step 5a, without Mont Blanc and
Taconnaz, worse by at most 0.05 C (F).

Usage:  forcing_runs_score.py [out_dir] [--runs today=ce25_pred_w0,F=ce25_F_w0,...]
"""
import fnmatch
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.validation import STEP5A_GLACIER_IDS

REPO = Path(__file__).resolve().parents[2]
RUNS = {'today': 'ce25_pred_w0', 'R': 'ce25_R_w0', 'F': 'ce25_F_w0', 'H': 'ce25_H_w0',
        'FH': 'ce25_FH_w0', 'P': 'ce25_P_w0'}
SUB = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')
CALI = Path('training_runs/monthly/CentralEurope/calibration/calibrate_m1_cID9_centraleurope_final_era5.dat')
CACHE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/cache/datahandler.pkl')
ARGS = [a for a in sys.argv[1:] if not a.startswith('--runs')]
for a in sys.argv[1:]:
    if a.startswith('--runs='):
        RUNS = dict(kv.split('=') for kv in a[len('--runs='):].split(','))
OUT = Path(ARGS[0] if ARGS else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/thermal_forcing/score')
DEEP, MB_FILES = 15.0, ('firnice_snowmay_*', 'firnice_snowbh_*', 'firnice_firnline_*')
MONT_BLANC, TACONNAZ = '00786', '00774'


def rmse(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else float('nan')


def run_dir(label):
    return REPO / 'data' / 'bayescal' / RUNS[label]


def files_equal(a, b, include='*', exclude=None):
    """(identical, differing names, only in one) of two runs' thermal output files."""
    keep = lambda n: fnmatch.fnmatch(n, include) and not (exclude and fnmatch.fnmatch(n, exclude))
    fa = {p.name for p in (run_dir(a) / SUB).iterdir() if keep(p.name)}
    fb = {p.name for p in (run_dir(b) / SUB).iterdir() if keep(p.name)}
    diff = [n for n in sorted(fa & fb)
            if (run_dir(a) / SUB / n).read_bytes() != (run_dir(b) / SUB / n).read_bytes()]
    return len(fa & fb) - len(diff), diff, sorted(fa ^ fb)


def toff_table():
    t = {}
    for line in (run_dir('today') / CALI).read_text().splitlines():
        p = line.split()
        if len(p) >= 12 and p[0].isdigit():
            t[p[0].zfill(5)] = float(p[11])
    return t


def points(labels, gl, toff):
    """One row per borehole point with the model of every run (finite in all of them)."""
    by_name = {g.glacier_name: g for g in gl}
    cg = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
    outs = {k: pickle.load(open(run_dir(k) / 'training_runs' / 'design_0000.output.pkl', 'rb')) for k in labels}
    rows = []
    for name in sorted(outs['today']):
        g = by_name.get(name)
        if g is None or any(outs[k].get(name) is None for k in labels):
            continue
        n = len(g.depths)
        vals = {k: np.asarray(outs[k][name][1], dtype=float)[:n] for k in labels}
        for i, d in enumerate(g.depths):
            if all(np.isfinite(vals[k][i]) for k in labels):
                rows.append(dict(entity=name, glacier=g.base_glacier_name, rgi=g.glacier_id, band=g.elevation_band,
                                 depth=float(d), obs=float(g.T_obs[i]), label=str(g.firn_label[i]),
                                 toff=toff.get(g.glacier_id, np.nan), test=g.glacier_id in STEP5A_GLACIER_IDS,
                                 cg_deep=(name, float(d)) in cg, **{k: float(vals[k][i]) for k in labels}))
    return pd.DataFrame(rows)


def perm_spearman(x, y, n=20000, seed=0):
    rho = stats.spearmanr(x, y).correlation
    rng = np.random.default_rng(seed)
    null = np.array([stats.spearmanr(x, rng.permutation(y)).correlation for _ in range(n)])
    return float(rho), float(np.mean(np.abs(null) >= abs(rho)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    labels = [k for k in RUNS if (run_dir(k) / 'training_runs' / 'design_0000.output.pkl').exists()]
    gl = pickle.load(open(CACHE, 'rb'))['glaciers']
    toff = toff_table()
    res, lines = {'runs': labels}, []
    say = lambda s='': (print(s), lines.append(s))

    # -- checks on files --------------------------------------------------------------------
    say('## Checks fixed before the runs\n')
    checks = []
    if 'R' in labels:
        same, diff, odd = files_equal('R', 'today', exclude='firnice_wetdepth_*')
        pa = pickle.load(open(run_dir('R') / 'training_runs/design_0000.output.pkl', 'rb'))
        pb = pickle.load(open(run_dir('today') / 'training_runs/design_0000.output.pkl', 'rb'))
        pk = all(k in pb and (pa[k] is None) == (pb[k] is None) and
                 (pa[k] is None or all(np.array_equal(np.asarray(x), np.asarray(y), equal_nan=True)
                                       for x, y in zip(pa[k], pb[k]))) for k in pa) and set(pa) == set(pb)
        checks.append(('R reproduces today', f'{same} files identical, {len(diff)} differ, {len(odd)} unmatched; '
                       f'matched values {"identical" if pk else "DIFFER"}', not diff and not odd and pk))
    for k in ('F', 'H', 'FH', 'P'):
        if k in labels:
            n_same = n_diff = 0
            for pat in MB_FILES:
                s, d, o = files_equal(k, 'today', include=pat)
                n_same += s
                n_diff += len(d) + len(o)
            checks.append((f'{k}: mass-balance files unchanged', f'{n_same} identical, {n_diff} differ', n_diff == 0))
    if 'F' in labels:
        zero = sorted(g for g, t in toff.items()
                      if t == 0.0 and (run_dir('F') / SUB / f'firnice_snowclim_{g}.dat').exists())
        n_same, bad = 0, []
        for g in zero:
            s, d, o = files_equal('F', 'today', include=f'*_{g}.dat', exclude='firnice_wetdepth_*')
            n_same += s
            bad += d + o
        checks.append((f'F: outlines with t_offset 0.00 ({", ".join(zero)}) byte-identical',
                       f'{n_same} identical, {len(bad)} differ', not bad))

    # -- points -------------------------------------------------------------------------------
    P = points(labels, gl, toff)
    P.to_csv(OUT / 'points.csv', index=False)
    fit = P[~P.cg_deep]
    deep = fit[fit.depth > DEEP]
    say(f'{P.entity.nunique()} entities, {P.glacier.nunique()} glaciers, {len(P)} points finite in all runs '
        f'({len(deep)} deeper than {DEEP:.0f} m outside Colle Gnifetti below 40 m)\n')

    # same site, two outlines
    cgs = P[P.rgi.isin(['01225', '01230']) & (P.band >= 4400) & (P.depth >= 10) & (P.depth <= 40)]
    same_site = {}
    for k in labels:
        b = {r: float((cgs[cgs.rgi == r][k] - cgs[cgs.rgi == r].obs).mean()) for r in ('01225', '01230')}
        same_site[k] = dict(bias_01225=b['01225'], bias_01230=b['01230'], diff=b['01230'] - b['01225'])
    for k in ('F', 'FH'):
        if k in labels:
            checks.append((f'{k}: Colle Gnifetti 01230 minus 01225 at 10-40 m below 1 C',
                           f"{same_site['today']['diff']:+.2f} -> {same_site[k]['diff']:+.2f} C", abs(same_site[k]['diff']) < 1))

    # t_offset against deep bias, per outline
    corr = {}
    for k in labels:
        per = deep.assign(e=deep[k] - deep.obs).groupby('rgi').agg(toff=('toff', 'first'), bias=('e', 'mean'))
        corr[k] = dict(zip(('rho', 'p'), perm_spearman(per.toff.to_numpy(), per.bias.to_numpy()))) | {'n': len(per)}
    for k in ('F', 'FH'):
        if k in labels:
            checks.append((f'{k}: deep bias no longer correlated with t_offset (permutation p > 0.05)',
                           f"rho {corr['today']['rho']:+.2f} (p {corr['today']['p']:.3f}) -> "
                           f"{corr[k]['rho']:+.2f} (p {corr[k]['p']:.3f}), {corr[k]['n']} outlines", corr[k]['p'] > 0.05))

    others = deep[~deep.test & ~deep.rgi.isin([MONT_BLANC, TACONNAZ])]
    if 'F' in labels:
        d_rmse = rmse(others.F - others.obs) - rmse(others.today - others.obs)
        checks.append(('F: deep RMSE outside Step 5a, without Mont Blanc and Taconnaz, worse by at most 0.05 C',
                       f'{rmse(others.today - others.obs):.3f} -> {rmse(others.F - others.obs):.3f} C ({d_rmse:+.3f})',
                       d_rmse <= 0.05))
    say('| check | result | pass |\n|---|---|---|')
    for name, txt, ok in checks:
        say(f'| {name} | {txt} | {"yes" if ok else "**no**"} |')
    res['checks'] = [dict(check=c, result=t, passed=bool(o)) for c, t, o in checks]
    res['same_site'], res['toff_correlation'] = same_site, corr

    # -- the table ----------------------------------------------------------------------------
    groups = [('all glaciers', np.ones(len(fit), bool)), ('test glaciers (Step 5a)', fit.test.to_numpy()),
              ('outside Step 5a', (~fit.test).to_numpy()),
              ('outside Step 5a, without Mont Blanc and Taconnaz',
               (~fit.test & ~fit.rgi.isin([MONT_BLANC, TACONNAZ])).to_numpy())]
    say('\n## RMSE / bias [C], points deeper than 15 m (and all points)\n')
    say('| glaciers | points | ' + ' | '.join(labels) + ' |\n|---|---|' + '---|' * len(labels))
    tab = []
    for gname, gm in groups:
        for dname, dm in (('> 15 m', (fit.depth > DEEP).to_numpy()), ('all', np.ones(len(fit), bool))):
            sub = fit[gm & dm]
            cells = [f'{rmse(sub[k] - sub.obs):.2f} / {np.mean(sub[k] - sub.obs):+.2f}' for k in labels]
            say(f'| {gname}, {dname} | {len(sub)} | ' + ' | '.join(cells) + ' |')
            tab.append(dict(glaciers=gname, depth=dname, n=len(sub),
                            **{k: dict(rmse=rmse(sub[k] - sub.obs), bias=float(np.mean(sub[k] - sub.obs))) for k in labels}))
    res['table'] = tab
    say('\n## Sites [C]: measured, then each run (points deeper than 15 m where available)\n')
    say('| site | points | measured | ' + ' | '.join(labels) + ' |\n|---|---|---|' + '---|' * len(labels))
    sites = []
    for name, m in (('Mont Blanc 00786', P.rgi == MONT_BLANC),
                    *[(f'Taconnaz {b} m', (P.rgi == TACONNAZ) & (P.band == b))
                      for b in sorted(P[P.rgi == TACONNAZ].band.unique())],
                    ('Tete Rousse 00777, > 15 m', (P.rgi == '00777') & (P.depth > DEEP)),
                    ('Colle Gnifetti 01225, 10-40 m', P.index.isin(cgs[cgs.rgi == '01225'].index)),
                    ('Colle Gnifetti 01230, 10-40 m', P.index.isin(cgs[cgs.rgi == '01230'].index)),
                    ('Colle Gnifetti below 40 m (apart)', P.cg_deep)):
        sub = P[m]
        if len(sub):
            say(f'| {name} | {len(sub)} | {sub.obs.mean():.2f} | ' + ' | '.join(f'{sub[k].mean():.2f}' for k in labels) + ' |')
            sites.append(dict(site=name, n=len(sub), obs=float(sub.obs.mean()), **{k: float(sub[k].mean()) for k in labels}))
    res['sites'] = sites
    say('\n## Firn and ice points, > 15 m: RMSE [C]\n')
    for lab in ('firn', 'ice', 'unknown'):
        sub = deep[deep.label == lab]
        say(f'- {lab} ({len(sub)}): ' + ', '.join(f'{k} {rmse(sub[k] - sub.obs):.2f}' for k in labels))
    if all(k in labels for k in ('F', 'H', 'FH')):
        inter = (deep.FH - deep.F) - (deep.H - deep.today)
        say(f'\ninteraction (FH - F) - (H - today) at the deep points: mean {inter.mean():+.3f}, '
            f'max |.| {inter.abs().max():.3f} C')
        res['interaction'] = dict(mean=float(inter.mean()), max_abs=float(inter.abs().max()))
    if 'F' in labels:
        per = deep.assign(d=deep.F - deep.today).groupby('rgi').agg(glacier=('glacier', 'first'), toff=('toff', 'first'),
                                                                     n=('d', 'size'), change=('d', 'mean'))
        say('\n## F minus today per outline, > 15 m [C], against t_offset\n')
        say('| outline | glacier | t_offset | points | change |\n|---|---|---|---|---|')
        for rgi, r in per.sort_values('toff').iterrows():
            say(f'| {rgi} | {r.glacier} | {r.toff:.2f} | {r.n} | {r.change:+.2f} |')
        res['f_change'] = per.reset_index().to_dict('records')
    (OUT / 'score.md').write_text('\n'.join(lines) + '\n')
    (OUT / 'score.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
