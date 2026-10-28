"""Free-air mass balance, Step 2: score the three options for the air below ERA5's ground at the stations
(development_history/plans/2026-10-28_freeair_step2_candidates.md, fixed before this script ran).

(a) free air; (b) free air above the ground, today's air below; (c) free air above the ground, below it
(1 - w) * free air + w * today's air with w = min(d / 400, 1), d = depth below the ground. Series of chapter 13 at
GloGEM's climate point of each station (free air: 'node' for the mountain stations, 'nearest' = the same point for
the valley stations; today's air: glogem_series.csv), each station's full overlap with its observations.

Writes <out>/step2_scores.csv and prints the table. Run: python -I freeair_step2_score.py [out_dir]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ST = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/freeair/station_checks/tables')
OBS_HIGH = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/airtemp_check/tables/'
                'bias_monthly_series.csv')
OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/mbfreeair/step2')
D_BLEND = 400.0
SEASONS = {'DJF': (12, 1, 2), 'MAM': (3, 4, 5), 'JJA': (6, 7, 8), 'SON': (9, 10, 11)}
NEAR_GROUND = ['Grand St-Bernard', 'Bocchetta delle Pisse', 'Aiguilles Rouges']


def series():
    cells = pd.read_csv(ST / 'station_cells.csv')
    cells = cells[cells.group.isin(['ch12', 'valley'])]
    fa = pd.read_csv(ST / 'freeair_series.csv')
    fa = fa[((fa.group == 'ch12') & (fa.variant == 'node')) | ((fa.group == 'valley') & (fa.variant == 'nearest'))]
    t2 = pd.read_csv(ST / 'glogem_series.csv')
    o1 = pd.read_csv(OBS_HIGH)[['station', 'year', 'month', 'T_obs']]
    o2 = pd.read_csv(ST / 'obs_monthly.csv')[['station', 'year', 'month', 'T_obs']]
    out = []
    for _, c in cells.iterrows():
        obs = (o1 if c.group == 'ch12' else o2)
        obs = obs[obs.station == c.station].groupby(['year', 'month'], as_index=False).T_obs.mean()
        a = fa[(fa.group == c.group) & (fa.station == c.station)][['year', 'month', 'T']].rename(columns={'T': 'FA'})
        b = t2[(t2.group == c.group) & (t2.station == c.station)][['year', 'month', 'T']].rename(columns={'T': 'T2'})
        m = obs.merge(a, on=['year', 'month']).merge(b, on=['year', 'month']).dropna()
        m['station'], m['z'], m['d'] = c.station, c.z, c.node_hclim - c.z
        out.append(m)
    return pd.concat(out, ignore_index=True)


def options(df):
    d = df.d.values
    w = np.clip(d / D_BLEND, 0.0, 1.0)
    df['a'] = df.FA
    df['b'] = np.where(d > 0, df.T2, df.FA)
    df['c'] = (1 - w) * df.FA + w * df.T2
    return df


def score(df):
    rows = []
    for st, g in df.groupby('station', sort=False):
        grp = 'valley' if g.d.iloc[0] > 0 else ('near the ground' if st in NEAR_GROUND else 'above 2800 m')
        for opt in 'abc':
            e = g[opt] - g.T_obs
            r = dict(station=st, group=grp, z=g.z.iloc[0], depth=round(g.d.iloc[0]), option=opt, n=len(g),
                     first=f'{g.year.min()}', last=f'{g.year.max()}', year=e.mean(),
                     rmse=float(np.sqrt((e ** 2).mean())))
            for s, ms in SEASONS.items():
                r[s] = e[g.month.isin(ms)].mean()
            rows.append(r)
    return pd.DataFrame(rows)


def main(out=OUT):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    S = score(options(series()))
    S.to_csv(out / 'step2_scores.csv', index=False, float_format='%.3f')
    cols = ['DJF', 'MAM', 'JJA', 'SON', 'year', 'rmse']
    print('model minus station [C]; options: (a) free air, (b) step at the ground, (c) gradual over 400 m')
    for grp in ['valley', 'near the ground', 'above 2800 m']:
        G = S[S.group == grp]
        print(f'\n{grp} ({G.station.nunique()} stations)')
        print(G.groupby('option')[cols].mean().round(2).to_string())
    print('\nvalley stations, option by option')
    print(S[S.group == 'valley'].pivot_table(index=['station', 'depth'], columns='option', values=['DJF', 'rmse'])
          .round(2).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:2])
