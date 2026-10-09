md(r'''
## 2. Data and rows per target

The module `target_analysis.py` builds every row from DataHandler's own pooling (CentralEurope,
"estimated" and "unassessed" profiles included) and reads the model at each reading's date
through the pipeline's matcher (`runner.match_model_to_obs`). Its cache was rebuilt from scratch
and checked in Stage 1; the checks are rerun here (a fraction of a second).
''')

code(r'''
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from IPython.display import display

REPO = Path.cwd().resolve().parent                     # the notebook runs in notebooks/
if not (REPO / 'scripts' / 'icetemp_analysis' / 'target_analysis.py').exists():
    REPO = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp')
sys.path.insert(0, str(REPO / 'scripts' / 'icetemp_analysis'))
import target_analysis as ta                                             # noqa: E402
from icetemp.calibration.thermal_structure import pressure_melting_point  # noqa: E402

STEP5A = ta.SCRATCH / 'step5a'          # Stage 1-3 outputs: cache/, stage2/, identifiability/, heldout/, checks/
pd.set_option('display.precision', 3, 'display.width', 220, 'display.max_columns', 30, 'display.max_colwidth', 80)

GOT = {}                                # numbers computed here; compared with the stage reports in section 8


def got(key, value):
    GOT[key] = float(value)


# figures: one colour per target in every figure, neutral ink otherwise
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#a3a29c', '#e6e5e1'
COL = {'T1': '#2a78d6', 'T2': '#eb6834', 'T3': '#1baf7a'}
plt.rcParams.update({'figure.dpi': 110, 'font.size': 9.5, 'axes.titlesize': 10, 'axes.titlelocation': 'left',
                     'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
                     'xtick.color': INK2, 'ytick.color': INK2, 'legend.frameon': False, 'lines.linewidth': 2,
                     'lines.solid_capstyle': 'round', 'text.color': INK})
''')

code(r'''
obs = ta.load_obs()
s2, r7 = ta.load_model('alps_s2'), ta.load_model('alps_r7')
common = ta.common_rows([s2, r7], obs)             # points with a model value in all 64 runs, and their T3 rows
SD = ta.prior_sd()                                 # prior sd of f_ref (U[0, 1]), f_ins, f_adv
theta = s2['theta']                                # the 32 runs: f_ref, f_ins, f_adv (identical in alps_r7)
keep = s2['design'].keep.to_numpy()                # f_ref >= 0.1: the 29 runs the smooth stand-ins use
ik = np.flatnonzero(keep)
print(f"{len(obs['glaciers'])} entities, {len(obs['readings'])} readings; same design in both ensembles: "
      f"{np.array_equal(s2['theta'], r7['theta'])}; prior sd {np.round(SD, 4)}; runs left out of the smooth "
      f"fits (f_ref < 0.1): {[int(i) for i in np.flatnonzero(~keep)]}")


def build_targets(mod, date_run=ta.BEST_RUN):
    """T1, T2 and T3 (the T3r rows) over the common rows, each with Sigma = D + S, and for T3 the
    month-to-month term K and the date spread V (V from run 13's monthly values, as in Stage 3)."""
    t1, t2, t3 = (ta.target(n, obs, mod, common) for n in ('T1', 'T2', 'T3r'))
    t1.Sigma = ta.target_error_cov(t1, obs)
    t2.Sigma = ta.target_error_cov(t2, obs)
    t3.Sigma = (ta.target_error_cov(t3, obs, seasonal=ta.seasonal_cov)
                + np.diag(ta.date_spread_var(t3.rows, obs, mod['Y3'][date_run])))
    t1.A2 = t1.A[:, common[0]]                     # T1 model = A2 @ model at the T2 rows
    return {'T1': t1, 'T2': t2, 'T3': t3}


TG = build_targets(s2)
''')

code(r'''
checks = ta.run_checks(obs, {'alps_s2': s2, 'alps_r7': r7})
c = checks
rm13 = c['T2_rmse_design13']
display(pd.DataFrame([
    ('observations: T2 = count-weighted mean of its T3 rows',
     f"{c['T2_obs_vs_T3']['max_diff_in_window']:.1e} °C on {c['T2_obs_vs_T3']['in_window']} points", c['T2_obs_vs_T3']['passes']),
    ('model: T2 = count-weighted mean of its T3 rows (alps_s2 / alps_r7)',
     f"{c['alps_s2_T2_model_vs_T3']['max_diff']:.1e} / {c['alps_r7_T2_model_vs_T3']['max_diff']:.1e} °C over "
     f"{c['alps_s2_T2_model_vs_T3']['values_compared']} values each",
     c['alps_s2_T2_model_vs_T3']['passes'] and c['alps_r7_T2_model_vs_T3']['passes']),
    ('T1 = weighted mean of T2, pooling matrix vs calibrator.py (observations / model)',
     f"{c['alps_s2_T1_vs_T2']['max_diff_obs']:.1e} / {c['alps_s2_T1_vs_T2']['max_diff_model']:.1e} °C",
     c['alps_s2_T1_vs_T2']['passes'] and c['alps_r7_T1_vs_T2']['passes']),
    ('T3: the matcher vs a copy of its time interpolation',
     f"{c['alps_s2_T3_matcher_vs_at_date']['max_diff']:.1e} °C", c['alps_s2_T3_matcher_vs_at_date']['passes']),
    ('design 13, T2 RMSE on the common points (alps_s2 / alps_r7)',
     f"{rm13['alps_s2']['rmse_joint']:.4f} / {rm13['alps_r7']['rmse_joint']:.4f} °C on {rm13['alps_s2']['joint_rows']} points",
     rm13['alps_s2']['passes'] and rm13['alps_r7']['passes']),
], columns=['check (Stage 1)', 'result', 'passes']).set_index('check (Stage 1)'))
print('all Stage 1 checks pass:', c['all_pass'])
got('s1.rmse_d13_s2', rm13['alps_s2']['rmse_joint']); got('s1.rmse_d13_r7', rm13['alps_r7']['rmse_joint'])
''')

code(r'''
t3lit = ta.target('T3', obs, s2, common)          # T3 as first specified (range-dated readings spread over months)
rows_tab = pd.DataFrame({
    'one row is': ['an entity: depth-weighted mean of its points', 'a point: 1 m depth bin',
                   'a point and month; a range-dated reading repeated in each of its months',
                   'a point and month; a range-dated reading kept whole'],
    'rows compared': [len(TG['T1'].y), len(TG['T2'].y), len(t3lit.y), len(TG['T3'].y)],
    'rows in the tables': [len(obs['T1']), len(obs['T2']), len(obs['T3']), len(obs['T3r'])]},
    index=pd.Index(['T1', 'T2', 'T3 as first specified', 'T3 (the T3r rows, used from here on)'], name='target'))
display(rows_tab)
n_rng = int((TG['T3'].rows.kind == 'range').sum())
n_rep = int(t3lit.rows.range_only.sum())
print(f'{n_rng} readings are dated only to a range of 4-25 months; as first specified they fill {n_rep} of the '
      f'{len(t3lit.y)} T3 rows, as T3r they are {n_rng} rows next to {int((TG["T3"].rows.kind == "month").sum())} month rows.')
for k, v in (('s1.rows_T1', len(TG['T1'].y)), ('s1.rows_T2', len(TG['T2'].y)), ('s1.rows_T3lit', len(t3lit.y)),
             ('s1.rows_T3r', len(TG['T3'].y)), ('s1.range_readings', n_rng), ('s1.range_rows_lit', n_rep)):
    got(k, v)
z_out = obs['T2'].depth.to_numpy(float)[~common[0]]
print(f'points never compared (no model value in some run): {len(z_out)}, at {z_out.min():.3f}–{z_out.max():.3f} m depth')
''')

code(r'''
r2 = TG['T2'].rows
rgi = obs['T2'].groupby('system_short').rgi.agg(lambda s: ' + '.join(sorted(set(s))))
per_sys = pd.concat([rgi.rename('RGI7'),
                     TG['T1'].rows.groupby('system_short').size().rename('T1 entities'),
                     r2.groupby('system_short').size().rename('T2 points'),
                     TG['T3'].rows.groupby('system_short').size().rename('T3 rows'),
                     t3lit.rows.groupby('system_short').size().rename('T3 as first specified')], axis=1)
per_sys.index.name = 'glacier system'
per_sys.loc['all'] = ['', *per_sys.iloc[:, 1:].sum()]
display(per_sys)
for s, row in per_sys.drop('all').iterrows():
    got(f's1.sys_T2_{s}', row['T2 points']); got(f's1.sys_T3r_{s}', row['T3 rows'])
bands = [b for _, _, b in ta.DEPTH_BINS]
zt = pd.crosstab(pd.Categorical(ta.depth_band(r2.depth), bands), r2.firn_label.to_numpy(), margins=True, margins_name='all')
zt.index.name, zt.columns.name = 'T2 points by depth', 'firn label'
display(zt)
for b in bands:
    got(f's1.T2_depth_{b}', (ta.depth_band(r2.depth) == b).sum())
for lab in ('firn', 'ice', 'unknown'):
    got(f's1.T2_label_{lab}', (r2.firn_label == lab).sum())
print(f"T2 points that average over several months / seasons / years: {(r2.n_months > 1).sum()} / "
      f"{(r2.n_seasons > 1).sum()} / {(r2.n_years > 1).sum()}; with range-dated readings: {r2.range_dated.sum()}")
got('s1.T2_mix_months', (r2.n_months > 1).sum()); got('s1.T2_mix_seasons', (r2.n_seasons > 1).sum())
''')

code(r'''
rbt = pd.DataFrame(checks['rmse_by_target']).T.astype({'rows': int, 'best_run': int})
rbt.index = pd.MultiIndex.from_tuples([tuple(i.split(' ')) for i in rbt.index], names=['ensemble', 'target'])
rbt = rbt.rename(columns={'design13': 'RMSE design 13 [°C]', 'best_run': 'best run', 'best': 'RMSE best run [°C]'})
display(rbt)
for (ens, tgt), row in rbt.iterrows():
    got(f's1.rmse_{ens}_{tgt}', row['RMSE design 13 [°C]'])
''')

md(r'''
**Reading.**
- **T3 as first specified repeats readings.** 168 readings dated only to a range of months fill
  1,914 of its 3,417 rows. Each repeat compares the *same* reading with a different month of the
  model, so the differences between those rows are the model's own seasonal changes, not misfits,
  and each such reading would weigh 4–25 times. **From here on, T3 means the T3r rows:** a
  range-dated reading is one row, compared with the model averaged over its months, and the
  model's spread over those months is added to its measurement error (V). The plan's 1,674 T3
  rows on `alps_r7` correspond to this count (1,671 here), not to T3 as first specified.
- **Monte Rosa dominates:** 422 of the 771 points (55 %) and 776 of the 1,671 T3 rows (46 %).
  Corvatsch (41 points, 301 T3 rows), Hohsaas and Hintereisferner gain most rows from T3: loggers
  read their points over many months.
- 5 points at 0.065–0.295 m depth are never compared: they lie above the model's first node.
- The RMSEs above are plain averages over each target's own rows (no weights), so they describe
  the rows, not the targets' quality; held-out RMSE (section 5) puts all targets on one ruler.
''')
