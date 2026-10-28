"""Figures of notebook 07 ("What the free air does to GloGEM") and of book chapter 15; each returns a figure.

Part 1, the air: air_at_stations, lapse_rates, air_change, area_below_ground. Data: the station tables of chapter 13
(week3/climate_products, week3/freeair/station_checks), the free-air table and Step 0's band table
(freeair_band_change.py).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
MB = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/mbfreeair')
TABLE = W3 / 'freeair' / 'freeair_CentralEurope.nc'
BANDS = MB / 'step0' / 'band_change.csv'
MONTHS = 'JFMAMJJASOND'
HYD_YEARS = (2000, 2019)

# the book's palette (development_history/figures/make_ch13.py)
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, GREEN, LIGHT = '#2a78d6', '#eb6834', '#2f9e6b', '#9cc3ef'
DIVERGING = LinearSegmentedColormap.from_list('cold_warm', [BLUE, '#ecebe8', ORANGE])
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def air_at_stations():
    """Monthly bias against the 7 stations above 2800 m, and seasonal bias at the 6 valley stations."""
    S = pd.read_csv(W3 / 'climate_products/tables/summary_above_2800m.csv')
    S = S[S.window == 'common all']
    C = pd.read_csv(W3 / 'freeair/station_checks/tables/s15_coldpool.csv')
    C = C[C.window == 'full overlap']
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.8), facecolor=SURFACE, sharey=True,
                                 gridspec_kw={'width_ratios': [2.2, 1]})
    for ax in (a1, a2):
        style(ax)
        ax.axhline(0, color=INK, lw=0.8)
    rows = [('ERA5 free air', 'node (each station)', BLUE, '-', 'free air'),
            ('GloGEM today', 'primary node', INK2, '--', 'today: 2 m air + lapse rate')]
    for prod, var, col, ls, lab in rows:
        r = S[(S['product'] == prod) & (S.variant == var)].iloc[0]
        y = [r[f'm{m:02d}'] for m in range(1, 13)]
        a1.plot(range(12), y, color=col, ls=ls, lw=2, marker='o', ms=5, label=f'{lab} (year {r.annual:+.1f} °C)')
    a1.set_xticks(range(12), list(MONTHS))
    a1.set_ylabel('model minus station [°C]', color=INK2)
    a1.set_title('7 stations above 2800 m (2010-2016)', color=INK, fontsize=10, loc='left')
    a1.legend(frameon=False, loc='lower right', labelcolor=INK)
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    for k, (prod, col, mk) in enumerate([('free air bilinear', BLUE, 'o'), ('GloGEM today', INK2, 's')]):
        sub = C[C['product'] == prod]
        y = [sub[s].mean() for s in seasons]
        a2.plot(np.arange(4) + (k - 0.5) * 0.18, y, ls='none', marker=mk, ms=8, color=col,
                markeredgecolor=SURFACE, markeredgewidth=2)
        a2.vlines(np.arange(4) + (k - 0.5) * 0.18, [sub[s].min() for s in seasons], [sub[s].max() for s in seasons],
                  color=col, lw=2)
    a2.set_xticks(range(4), seasons)
    a2.set_title('6 valley stations, 400-1140 m\nbelow ERA5\'s ground (mean, range)', color=INK, fontsize=10,
                 loc='left')
    fig.tight_layout()
    return fig


def _node_profiles():
    """Free-air climatology per climate point used by a glacier: (nodes, 12, heights), and the heights."""
    df = pd.read_csv(BANDS, usecols=['node_lon', 'node_lat'])
    nodes = df.drop_duplicates().values
    tab = xr.open_dataset(TABLE)
    hyd = tab.year.values.astype(int) + (tab.month.values.astype(int) >= 10)
    win = (hyd >= HYD_YEARS[0]) & (hyd <= HYD_YEARS[1])
    mon = tab.month.values.astype(int)
    prof = []
    for lon, lat in nodes:
        fa = tab.tfa.sel(longitude=lon, latitude=lat, method='nearest').values.astype(np.float64)
        prof.append([fa[win & (mon == m)].mean(axis=0) for m in range(1, 13)])
    return nodes, np.array(prof), tab.height.values.astype(float)


def lapse_rates():
    """The free air's own lapse rate by month in three height layers, against GloGEM's fixed dtdz (medians over
    the climate points of all glaciers)."""
    import freeair_band_change as fbc
    nodes, prof, hgt = _node_profiles()
    clim = fbc.read_clim()
    il = [int(np.argmin(np.abs(clim['rlon'] - lon))) for lon, _ in nodes]
    ia = [int(np.argmin(np.abs(clim['rlat'] - lat))) for _, lat in nodes]
    dtdz = np.array([clim['rtg'][:, i, j] for i, j in zip(il, ia)]) * 100      # K per 100 m
    fig, ax = plt.subplots(figsize=(7.2, 4.4), facecolor=SURFACE)
    style(ax)
    for (lo, hi), col in zip([(1500, 2500), (2500, 3500), (3500, 4500)], [LIGHT, BLUE, '#1a4f8f']):
        k0, k1 = int(np.argmin(np.abs(hgt - lo))), int(np.argmin(np.abs(hgt - hi)))
        g = (prof[:, :, k1] - prof[:, :, k0]) / (hi - lo) * 100
        ax.plot(range(12), np.median(g, axis=0), color=col, lw=2, marker='o', ms=5, label=f'free air {lo}-{hi} m')
    ax.plot(range(12), np.median(dtdz, axis=0), color=ORANGE, lw=2, ls='--', marker='s', ms=5,
            label='GloGEM today: dtdz (fixed per month)')
    ax.set_xticks(range(12), list(MONTHS))
    ax.set_ylabel('lapse rate [°C per 100 m]', color=INK2)
    ax.set_title(f'Median over the {len(nodes)} climate points of the Alpine glaciers, '
                 f'{HYD_YEARS[0]}-{HYD_YEARS[1]}', color=INK, fontsize=10, loc='left')
    ax.legend(frameon=False, labelcolor=INK, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
    fig.tight_layout()
    return fig


def air_change(step=250):
    """Free air minus today's band air (without T_off), area-weighted, by band height and month."""
    df = pd.read_csv(BANDS)
    d = df[[f'tfa_{m:02d}' for m in range(1, 13)]].values - df[[f't2m_{m:02d}' for m in range(1, 13)]].values
    a = df.area.values
    edges = np.arange(1250, 4750 + step, step)
    rows, labs, areas = [], [], []
    for lo in edges[:-1]:
        s = (df.z.values >= lo) & (df.z.values < lo + step)
        if a[s].sum() < 0.5:
            continue
        rows.append((d[s] * a[s, None]).sum(axis=0) / a[s].sum()); labs.append(f'{lo}-{lo + step} m')
        areas.append(a[s].sum())
    M = np.array(rows)
    fig, ax = plt.subplots(figsize=(8.4, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    im = ax.imshow(M, cmap=DIVERGING, norm=TwoSlopeNorm(0, vmin=-np.abs(M).max(), vmax=np.abs(M).max()),
                   aspect='auto', origin='lower')
    for i in range(M.shape[0]):
        for j in range(12):
            ax.text(j, i, f'{M[i, j]:+.1f}', ha='center', va='center', fontsize=7.5, color=INK)
    ax.set_xticks(range(12), list(MONTHS))
    ax.set_yticks(range(len(labs)), [f'{l}  ({ar:.0f} km²)' for l, ar in zip(labs, areas)])
    ax.tick_params(colors=INK2, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label('free air minus today\'s air [°C]', color=INK2)
    cb.outline.set_visible(False)
    ax.set_title(f'All {df.id.nunique()} Alpine glaciers, {HYD_YEARS[0]}-{HYD_YEARS[1]}, area-weighted',
                 color=INK, fontsize=10, loc='left')
    fig.tight_layout()
    return fig


def area_below_ground():
    """Glacier area below ERA5's ground, by depth below it, with the depths of the valley stations."""
    df = pd.read_csv(BANDS, dtype={'id': str})
    dz = df.hclim.values - df.z.values
    bins = [0, 100, 200, 300, 500, 800, 1200]
    ar = [df.area.values[(dz > lo) & (dz <= hi)].sum() for lo, hi in zip(bins[:-1], bins[1:])]
    C = pd.read_csv(W3 / 'freeair/station_checks/tables/s15_coldpool.csv').drop_duplicates('station')
    fig, ax = plt.subplots(figsize=(7.2, 3.4), facecolor=SURFACE)
    style(ax)
    y = np.arange(len(ar))
    ax.barh(y, ar, height=0.7, color=BLUE)
    for yi, v in zip(y, ar):
        ax.text(v + 0.2, yi, f'{v:.1f} km²', va='center', color=INK, fontsize=9)
    ax.set_yticks(y, [f'{lo}-{hi} m' for lo, hi in zip(bins[:-1], bins[1:])])
    ax.invert_yaxis()
    ax.set_xlabel('glacier area [km²]', color=INK2)
    ax.set_ylabel('depth below ERA5\'s ground', color=INK2)
    C['k'] = np.searchsorted(bins, C.below_era5_ground_m.values) - 1
    for k, g in C.groupby('k'):
        ax.text(max(ar) * 0.55, k, 'stations: ' + ', '.join(g.station), va='center', fontsize=8, color=INK2)
    tot = df.area.sum()
    ax.set_title(f'{sum(ar):.1f} km² of {tot:.0f} km² ({100 * sum(ar) / tot:.1f} %) lie below ERA5\'s ground; '
                 'valley stations in grey', color=INK, fontsize=10, loc='left')
    fig.tight_layout()
    return fig
