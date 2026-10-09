"""Stage 3a tables: CSV files and tables.md from results.pkl (s3a_main.py)."""
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
R = pickle.load(open(OUT / 'results.pkl', 'rb'))
TNAME = {'T1': 'T1', 'T2': 'T2', 'T3r': 'T3 (=T3r)', 'T3lit': 'T3 literal*'}
ROWS = {}


def f(x, nd=2):
    try:
        if x is None or not np.isfinite(float(x)):
            return '–'
    except (TypeError, ValueError):
        return str(x)
    return f'{float(x):.{nd}f}'


def md(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|']
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    return '\n'.join(out)


def rng(v, nd=2):
    v = np.asarray(v, float)
    return f'{v.min():.{nd}f}–{v.max():.{nd}f}' if len(v) else '–'


def weak(r):
    w = [float(x) for x in str(r.weak_dir).split()]
    names = ['f_ref', 'f_ins', 'f_adv'][:len(w)]
    return ' '.join(f'{x:+.2f} {n}' for x, n in zip(w, names)) + f' (sd {r.weak_sd:.2f})'


def main_table(ens, mode, targets, variant='primary'):
    main, ver = R['main'], R['verdicts']
    out = []
    for t in targets:
        for meth in ('surface', 'emulator'):
            m = main[(main.ens == ens) & (main.target == t) & (main['mode'] == mode) & (main.method == meth)
                     & (main.variant == variant)]
            if m.empty:
                continue
            b = m[m.point == 'best point'].iloc[0]
            runs = m[m.point.str.startswith('run ')]
            v = ver[(ver.ens == ens) & (ver.target == t) & (ver['mode'] == mode) & (ver.method == meth)
                    & (ver.variant == variant)].iloc[0]
            row = {'target': TNAME[t], 'method': meth,
                   'best point': f'{b.f_ref:.2f}, {b.f_ins:.2f}' + (f', {b.f_adv:.2f}' if mode == '3p' else ''),
                   'corr': f(b.corr_ref_ins), 'spread f_ref': f(b.spread_f_ref), 'spread f_ins': f(b.spread_f_ins)}
            if mode == '3p':
                row['spread f_adv'] = f(b.spread_f_adv)
            row['cond'] = f(b.cond_post, 1)
            row['weakest direction'] = weak(b)
            row['5 best runs'] = v.runs
            row['runs: corr'] = rng(runs.corr_ref_ins)
            row['runs: spread f_ins'] = rng(runs.spread_f_ins)
            if mode == '3p':
                row['runs: spread f_adv'] = rng(runs.spread_f_adv)
            row['runs passing'] = f'{v.runs_pass}/{v.n_runs}'
            row['told apart'] = 'PASS' if v.told_apart else 'fail'
            if mode == '3p':
                row['f_adv rule'] = 'PASS' if v.adv_rule else 'fail'
            if meth == 'emulator' and 'E_mean_var' in b and np.isfinite(b.E_mean_var):
                row['E at best point (mean var / top eig, C^2)'] = f'{b.E_mean_var:.3f} / {b.E_max_eig:.2f}'
            out.append(row)
    return pd.DataFrame(out).fillna('')


def ref_table():
    main = R['main']
    m = main[(main.point.str.startswith('reference')) & (main.variant == 'primary') & (main.target != 'T3lit')]
    out = []
    for _, r in m.iterrows():
        out.append({'ensemble': r.ens, 'target': TNAME[r.target], 'mode': r['mode'], 'method': r.method,
                    'corr': f(r.corr_ref_ins), 'spread f_ref': f(r.spread_f_ref), 'spread f_ins': f(r.spread_f_ins),
                    'spread f_adv': f(r.spread_f_adv) if r['mode'] == '3p' else '', 'cond': f(r.cond_post, 1)})
    return pd.DataFrame(out)


def sens_table(mode):
    sens, ver, main = R['sens'], R['verdicts'], R['main']
    out = []
    for meth in ('surface', 'emulator'):
        row = {'variant': 'primary', 'method': meth}
        for t in ('T1', 'T2', 'T3r'):
            b = main[(main.ens == 'alps_s2') & (main.target == t) & (main['mode'] == mode) & (main.method == meth)
                     & (main.variant == 'primary') & (main.point == 'best point')].iloc[0]
            v = ver[(ver.ens == 'alps_s2') & (ver.target == t) & (ver['mode'] == mode) & (ver.method == meth)
                    & (ver.variant == 'primary')].iloc[0]
            txt = f'{f(b.corr_ref_ins)} / {f(b.spread_f_ins)}' + (f' / {f(b.spread_f_adv)}' if mode == '3p' else '')
            row[TNAME[t]] = txt + f' ({v.runs_pass}/5, {"PASS" if v.told_apart else "fail"}' + \
                (f', f_adv {"PASS" if v.adv_rule else "fail"})' if mode == '3p' else ')')
        out.append(row)
    for vname in dict.fromkeys(sens.variant):
        for meth in ('surface', 'emulator'):
            row = {'variant': vname, 'method': meth}
            ok = False
            for t in ('T1', 'T2', 'T3r'):
                m = sens[(sens.variant == vname) & (sens.method == meth) & (sens.target == t) & (sens['mode'] == mode)]
                if m.empty:
                    continue
                ok = True
                b = m[m.point == 'best point'].iloc[0]
                v = ver[(ver.variant == vname) & (ver.method == meth) & (ver.target == t) & (ver['mode'] == mode)].iloc[0]
                txt = f'{f(b.corr_ref_ins)} / {f(b.spread_f_ins)}' + (f' / {f(b.spread_f_adv)}' if mode == '3p' else '')
                row[TNAME[t]] = txt + f' ({v.runs_pass}/5, {"PASS" if v.told_apart else "fail"}' + \
                    (f', f_adv {"PASS" if v.adv_rule else "fail"})' if mode == '3p' else ')')
            if ok:
                out.append(row)
    return pd.DataFrame(out)


def info_table(target, grouping, method='surface', mode='2p'):
    i = R['info']
    m = i[(i.target == target) & (i.grouping == grouping) & (i.method == method) & (i['mode'] == mode)]
    out = []
    for _, r in m.iterrows():
        row = {'group': r.group, 'rows': r.rows, 'share f_ref': f'{100 * r.share_f_ref:.0f} %',
               'share f_ins': f'{100 * r.share_f_ins:.0f} %', 'share along weakest direction': f'{100 * r.share_weak:.0f} %'}
        if mode == '3p':
            row['share f_adv'] = f'{100 * r.share_f_adv:.0f} %'
        row['alone: corr / spread f_ins'] = f'{f(r.alone_corr_ref_ins)} / {f(r.alone_spread_f_ins)}'
        row['without: corr / spread f_ins'] = (f'{f(r.without_corr_ref_ins)} / {f(r.without_spread_f_ins)}'
                                               if np.isfinite(r.get('without_corr_ref_ins', np.nan)) else '–')
        if mode == '3p':
            row['alone / without: spread f_adv'] = f'{f(r.alone_spread_f_adv)} / {f(r.without_spread_f_adv)}'
        out.append(row)
    return pd.DataFrame(out)


def quality_table():
    q = R['quality']
    out = []
    for _, r in q.iterrows():
        out.append({'target': TNAME[r.target], 'fitted on': r.fitted_on, 'depth': r.depth, 'rows': r.rows,
                    'flat rows (design sd < 0.02 C)': r.flat_rows, 'R2 median': f(r.r2_median, 3), 'R2 p10': f(r.r2_p10, 3),
                    'R2 p1': f(r.r2_p1, 3), 'R2 min': f(r.r2_min, 3), 'share R2 < 0.9': f'{100 * r.share_r2_below_0_9:.1f} %',
                    'R2 pooled': f(r.r2_pooled, 4),
                    'LOO RMSE surface (median / p90 / all), C': (f'{f(r.surface_loo_rmse_median, 3)} / {f(r.surface_loo_rmse_p90, 3)} / '
                                                                f'{f(r.surface_loo_rmse_all, 3)}' if np.isfinite(r.get('surface_loo_rmse_median', np.nan)) else '–'),
                    'LOO RMSE emulator (median / p90 / all), C': (f'{f(r.emulator_loo_rmse_median, 3)} / {f(r.emulator_loo_rmse_p90, 3)} / '
                                                                 f'{f(r.emulator_loo_rmse_all, 3)}' if np.isfinite(r.get('emulator_loo_rmse_median', np.nan)) else '–')})
    return pd.DataFrame(out)


def write():
    for k, name in (('main', 'fisher_all_points.csv'), ('verdicts', 'verdicts.csv'), ('sens', 'sensitivities.csv'),
                    ('info', 'information_sources.csv'), ('quality', 'surface_quality.csv'), ('jagree', 'j_agreement.csv'),
                    ('svd', 'svd_truncation.csv'), ('bestpts', 'best_points.csv'), ('runs', 'runs_ranked.csv'),
                    ('grid', 'grid_posterior_check.csv'), ('emu_meta', 'emulator_meta.csv'), ('emu_own', 'emulator_own_map.csv'),
                    ('decomp', 'emulator_objective_decomposition.csv')):
        R[k].to_csv(OUT / name, index=False, float_format='%.4f')
    L = ['# Stage 3a tables (generated by s3b_report.py from results.pkl)\n',
         'corr = posterior correlation f_ref-f_ins; spread = posterior sd / prior sd (1 = nothing learned); '
         'cond = largest / smallest eigenvalue of the posterior precision; weakest direction = unit vector in '
         'prior-sd units with its posterior sd. TOLD APART: spread f_ins <= 0.5 and |corr| <= 0.8 at the best point '
         'and at >= 4 of the 5 best runs. * literal T3 fails the honesty check (Stage 2) and is shown only for reference.\n']
    for ens in ('alps_s2', 'alps_r7'):
        tg = ('T1', 'T2', 'T3r') + (('T3lit',) if ens == 'alps_s2' else ())
        for mode in ('2p', '3p'):
            L.append(f'\n## {ens}: {"f_adv fixed at 1.0 (2 parameters)" if mode == "2p" else "f_adv free (3 parameters)"}\n')
            L.append(md(main_table(ens, mode, tg)))
    L.append('\n## Reference point: run 13 (2p at f_adv = 1.0, 3p at its own f_adv 1.055)\n')
    L.append(md(ref_table()))
    for mode in ('2p', '3p'):
        L.append(f'\n## Sensitivities, alps_s2, {mode}: corr / spread f_ins{" / spread f_adv" if mode == "3p" else ""} '
                 'at the best point (runs passing, verdict)\n')
        L.append(md(sens_table(mode)))
    for t, groupings in (('T1', ('glacier system', 'firn label (entity majority)')),
                         ('T2', ('glacier system', 'depth band', 'firn label', 'season', 'system x depth band')),
                         ('T3r', ('glacier system', 'depth band', 'firn label', 'season', 'depth band x season',
                                  'system x depth band'))):
        for g in groupings:
            for meth in ('surface', 'emulator'):
                L.append(f'\n## Where the information comes from: {TNAME[t]}, by {g} ({meth}, 2p, at the best point)\n')
                L.append(md(info_table(t, g, meth, '2p')))
        for g in ('glacier system', 'depth band', 'firn label') if t != 'T1' else ('glacier system',):
            L.append(f'\n## f_adv information: {TNAME[t]}, by {g} (surface, 3p, at the best point)\n')
            L.append(md(info_table(t, g, 'surface', '3p')))
    L.append('\n## Response-surface quality (quadratic, 29 runs) and emulator leave-one-out accuracy, alps_s2\n')
    L.append(md(quality_table()))
    L.append('\n## Sensitivities of the two methods compared at the best point (2p), alps_s2\n')
    L.append(md(R['jagree'].round(3).fillna('')))
    L.append('\n## SVD truncation: surface J projected on the emulator basis, alps_s2 (2p, best point)\n')
    L.append(md(R['svd'].round(3).fillna('')))
    L.append('\n## Best points (MAP inside the design-coverage gate; nn = distance to the nearest run, gate 1.26)\n')
    L.append(md(R['bestpts'].round(3)))
    L.append('\n## The emulator at its own MAP (diagnostic, not used in the verdict)\n')
    eo = R['emu_own'][['ens', 'target', 'mode', 'point', 'f_ref', 'f_ins', 'f_adv', 'corr_ref_ins', 'spread_f_ref',
                       'spread_f_ins', 'spread_f_adv', 'E_mean_var', 'E_max_eig']].round(3).fillna('')
    L.append(md(eo))
    L.append('\n## Emulator objective split at three points (-log posterior = half chi2 + half logdet - log prior)\n')
    L.append(md(R['decomp'].round(2)))
    L.append('\n## Grid posterior check (2p, exact posterior on a 0.01 grid inside the gate; spreads in prior sd)\n')
    L.append(md(R['grid'].round(3)))
    L.append('\n## Emulator meta (p components, truncation variance, LOO z std, inflation curve)\n')
    L.append(md(R['emu_meta']))
    L.append('\n## Runs ranked by the GLS misfit (top 8 per target)\n')
    rr = R['runs']
    L.append(md(rr[rr['rank'] <= 8].round(3)))
    (OUT / 'tables.md').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print(f'wrote {OUT / "tables.md"} and the CSV files')


if __name__ == '__main__':
    write()
