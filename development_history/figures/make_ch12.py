"""Figures for chapter 12 (sun and shade, burial, air: three tests and the orientation build).
Run from this folder: python make_ch12.py [preview_dir]

ch12_shade_slope: the shade error inside sites (within-site slope of the band bias on summer sunshine, OLS with site
  fixed effects; 90 % interval from the wild cluster bootstrap; all without the three entities Stage 1 left out of
  the within-site measures) for FRHS, Stage 1 (C, W, S), with the measured
  burial (B, B+C) and the build (three parts), against the rule-B line (half of the FRHS slope).
ch12_burial: the thermal model's burial rate (b_ref) against the measured long-term accumulation at the 41 borehole
  bands with a measurement (week3/burial_tables/burial.csv).
ch12_air_bias: the model's air temperature minus the station's, by month, at the high stations of Monte Rosa, Mont
  Blanc and the Bernese Alps / Engadin (week3/airtemp_check/tables/bias_by_station.csv, full overlap).
ch12_parts: the part factors F (central, dT 0) by elevation on the Gorner complex (01225) and Fieschergletscher
  (02611), with each part's share of the band (week3/orient_parts_v2/parts.csv).
ch12_results: the error on left-out massifs overall, in the zone and above 3800 m for FRHS, Stage 1, the measured
  burial, burial with orientation, and the build.
ch12_eb_day: the energy balance behind F on an average clear July day (Jungfraujoch days with at least 80 % of the
  possible sunshine, 1981-2020) at the example band (Aletschgletscher 02596, 3505 m), for flat open ground and
  30-deg planes facing north, east, south and west without a terrain horizon: the absorbed sun of each plane against
  the sun needed before melt starts (-psi), and the melt per hour. The model's own code (part_factors.Sun, eb_jjas).
ch12_eb_factor: F of the same idealised planes (June-September 1981-2020, the full calculation) against slope at
  dT 0, and against dT at 30 deg; beside it the band's real parts with three and with five parts.
ch12_pixels: every pixel's F at dT 0 between 3300 and 3700 m on the 23 calibration glaciers against its aspect, by
  slope class, with the part boundaries of three and five parts (week3/orient_pixels, part_factors.py --pixels).
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
ROOT = Path(__file__).resolve().parents[2]
ORI = ROOT / 'scripts/icetemp_analysis/orientation'
MTH = ROOT / 'data/bayescal/mth_FRHS_r0.85_i1.017_a1.00_mth_w0'
EX_GID, EX_BAND, EX_NAME = '02596', 3505.0, 'Aletschgletscher (02596), 3505 m'
EX_LAT = 46.50   # the glacier's centre; sun positions as in part_factors.py
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, NEUTRAL, LIGHT, GREEN = '#2a78d6', '#eb6834', '#9a9893', '#9cc3ef', '#2f9e6b'
PURPLE = '#8e5bd0'
PART_COL = {'G': NEUTRAL, 'N': BLUE, 'S': ORANGE, 'E': GREEN, 'W': PURPLE}
PLANES = (('flat', 0.0, 0.0, 'G'), ('north', 30.0, 0.0, 'N'), ('east', 30.0, 90.0, 'E'), ('south', 30.0, 180.0, 'S'),
          ('west', 30.0, 270.0, 'W'))


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg', facecolor=SURFACE)
    if PREVIEW:
        fig.savefig(PREVIEW / f'{name}.png', dpi=100, facecolor=SURFACE)
    plt.close(fig)


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def load(path):
    return json.load(open(path))['variants']


def shade_slope():
    s1, rc = load(W3 / 'orient_stage1/orient_stage1.json'), load(W3 / 'orient_recheck/orient_stage1.json')
    cp = load(W3 / 'orient_parts_score_v2/parts_score.json')
    rows = [('FRHS', s1['FRHS']['within']), ('Stage 1, central', s1['C']['within']),
            ('Stage 1, weaker', s1['W']['within']), ('Stage 1, stronger', s1['S']['within']),
            ('measured burial', rc['FRHS']['within']), ('burial + orientation', rc['C']['within']),
            ('the build (three parts)', cp['CP']['within (Stage-1 exclusions)'])]
    fig, ax = plt.subplots(figsize=(7.2, 3.8), facecolor=SURFACE)
    style(ax)
    y = np.arange(len(rows))[::-1]
    for yi, (lab, w) in zip(y, rows):
        lo, hi = w['slope_90']
        col = NEUTRAL if lab in ('FRHS', 'measured burial') else BLUE
        ax.plot([lo, hi], [yi, yi], color=col, lw=2, alpha=0.5)
        ax.plot(w['slope'], yi, 'o', color=col, ms=7)
    half = 0.5 * s1['FRHS']['within']['slope']
    ax.axvline(half, color=ORANGE, lw=1.2, ls='--')
    ax.axvline(0, color=INK2, lw=0.8)
    ax.text(half, len(rows) - 0.4, ' rule B: at least half of FRHS', color=ORANGE, fontsize=9, va='bottom')
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel('shade error inside sites (°C per +0.1 of summer sunshine)')
    ax.set_ylim(-0.7, len(rows) - 0.1)
    fig.tight_layout()
    save(fig, 'ch12_shade_slope')


def burial():
    B = pd.read_csv(W3 / 'burial_tables/burial.csv', dtype={'rgi': str})
    massif = {'00773': 'Mont Blanc', '00774': 'Mont Blanc', '01195': 'Monte Rosa', '01205': 'Monte Rosa',
              '01225': 'Monte Rosa', '01230': 'Monte Rosa'}
    B['massif'] = B.rgi.map(massif).fillna('other')
    fig, ax = plt.subplots(figsize=(6.4, 5.0), facecolor=SURFACE)
    style(ax)
    for m, col in (('Mont Blanc', ORANGE), ('Monte Rosa', BLUE), ('other', NEUTRAL)):
        d = B[B.massif == m]
        ax.errorbar(d.acc_central, d.model_bref, xerr=[d.acc_central - d.acc_low, d.acc_high - d.acc_central],
                    fmt='o', color=col, ms=5, lw=0.8, alpha=0.85, label=m)
    ax.plot([0, 6.5], [0, 6.5], color=INK2, lw=0.8)
    ax.text(3.0, 3.25, 'model = measured', color=INK2, fontsize=9, rotation=28, ha='center')
    ax.set_xlim(0, 3.6)
    ax.set_ylim(0, 6.5)
    ax.set_xlabel('measured long-term accumulation (m w.e. per year)')
    ax.set_ylabel("the thermal model's burial rate (m w.e. per year)")
    ax.legend(frameon=False, loc='upper right')
    fig.tight_layout()
    save(fig, 'ch12_burial')


def air_bias():
    A = pd.read_csv(W3 / 'airtemp_check/tables/bias_by_station.csv', dtype={'node': str})
    A = A[A.window == 'full overlap'].drop_duplicates('station')
    keep = {'Capanna Margherita': BLUE, 'Plateau Rosa': BLUE, 'Gornergrat': BLUE, 'Aiguille du Midi': ORANGE,
            'Jungfraujoch': NEUTRAL, 'Corvatsch': NEUTRAL}
    ls = {'Capanna Margherita': '-', 'Plateau Rosa': '--', 'Gornergrat': ':', 'Aiguille du Midi': '-',
          'Jungfraujoch': '-', 'Corvatsch': '--'}
    fig, ax = plt.subplots(figsize=(7.0, 3.8), facecolor=SURFACE)
    style(ax)
    months = np.arange(1, 13)
    for st, col in keep.items():
        r = A[A.station == st]
        if r.empty:
            continue
        r = r.iloc[0]
        ax.plot(months, [r[f'm{m:02d}'] for m in months], color=col, ls=ls[st], lw=1.8,
                label=f'{st} ({int(r.z)} m): year {r.annual:+.1f} °C')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(months)
    ax.set_xticklabels(list('JFMAMJJASOND'))
    ax.set_ylabel('model air minus measured (°C)')
    ax.legend(frameon=False, fontsize=8, loc='lower center', ncol=2)
    fig.tight_layout()
    save(fig, 'ch12_air_bias')


def parts():
    R = pd.read_csv(W3 / 'orient_parts_v2/parts.csv', dtype={'rgi': str})
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.4), sharey=False, facecolor=SURFACE)
    for ax, (gid, name) in zip(axes, (('01225', 'Gorner complex (01225)'), ('02611', 'Fieschergletscher (02611)'))):
        style(ax)
        d = R[R.rgi == gid]
        for p, lab in (('G', 'gentle'), ('N', 'north-facing'), ('S', 'south-facing')):
            q = d[(d.part == p) & (d.flag.isna())].sort_values('band')
            ax.scatter(q['F+0'], q.band, s=8 + 60 * q.share, color=PART_COL[p], alpha=0.75, label=lab)
        ax.axvline(1, color=INK2, lw=0.8)
        ax.set_title(name, fontsize=10, color=INK)
        ax.set_xlabel('orientation factor F (dot size: share of the band)')
    axes[0].set_ylabel('band elevation (m)')
    axes[0].legend(frameon=False, loc='lower left', fontsize=9)
    fig.tight_layout()
    save(fig, 'ch12_parts')


def results():
    s1, rc = load(W3 / 'orient_stage1/orient_stage1.json'), load(W3 / 'orient_recheck/orient_stage1.json')
    cp = load(W3 / 'orient_parts_score_v2/parts_score.json')
    rows = [('FRHS', s1['FRHS']['lomo']), ('Stage 1', s1['C']['lomo']), ('measured burial', rc['FRHS']['lomo']),
            ('burial + orientation', rc['C']['lomo']), ('the build', cp['CP']['lomo'])]
    keys = (('M1', 'all bands'), ('zone_M1', 'zone 3000-3800 m'), ('above_M1', 'above 3800 m'))
    fig, ax = plt.subplots(figsize=(7.4, 3.6), facecolor=SURFACE)
    style(ax)
    x = np.arange(len(keys))
    w = 0.16
    cols = [NEUTRAL, LIGHT, '#c9c7c2', BLUE, GREEN]
    for i, ((lab, L), col) in enumerate(zip(rows, cols)):
        vals = [L[k] for k, _ in keys]
        ax.bar(x + (i - 2) * w, vals, w, color=col, label=lab)
        for xi, v in zip(x, vals):
            ax.text(xi + (i - 2) * w, v + 0.03, f'{v:.2f}', ha='center', fontsize=7, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([k[1] for k in keys])
    ax.set_ylabel('error on left-out massifs (°C)')
    ax.legend(frameon=False, fontsize=8, ncol=3, loc='upper left')
    ax.set_ylim(0, 3.6)
    fig.tight_layout()
    save(fig, 'ch12_results')


class Example:
    """The example band's weather and energy balance, with the model's own code (part_factors.py)."""

    def __init__(self):
        sys.path[:0] = [str(ORI), str(ORI.parent)]
        import exposure_check as ec
        import orientation_eb as oe
        import part_factors as pf
        import write_orient_tables as wot
        self.ec, self.oe, self.pf = ec, oe, pf
        inp = pf.band_inputs(MTH, EX_GID, np.array([EX_BAND]))[EX_BAND]
        w = oe.read_station(wot.INPUTS['jungfraujoch'])
        sig, dif, _ = oe.sky_weights(w)
        from scipy.optimize import brentq
        R = pd.read_csv(W3 / 'orient_parts_v2/parts.csv', dtype={'rgi': str})
        p0_table = float(R[(R.rgi == EX_GID) & (R.band == EX_BAND)].psi0.iloc[0])
        self.dd = float(inp['dd_jjas'])
        self.lat = EX_LAT
        self.Th = oe.Forcing(inp['Tm'], w, sig, dif).Th[:, pf.JJAS_DAYS, :]
        self.sig = sig.reshape(oe.NY, 365)[:, pf.JJAS_DAYS, None]
        self.dif = dif.reshape(oe.NY, 365)[:, pf.JJAS_DAYS, None]
        self.sun = pf.Sun(self.lat)
        self.beam, fdif, fdir = self.sun.beam(EX_BAND)
        self.fdif = pf.hour_means(fdif)
        self.S_flat = self.sig * pf.hour_means(fdir)[None] + self.dif * self.fdif[None]
        target = max(self.dd, pf.DMIN)
        self.p0 = brentq(lambda x: pf.eb_jjas(self.Th, self.S_flat, x) - target, -900, 400, xtol=0.02)
        self.flat = {dT: pf.eb_jjas(self.Th + dT, self.S_flat, self.p0) for dT in pf.DTS}
        print(f'  example band: June-September degree-day melt {self.dd:.0f} mm, psi0 {self.p0:.2f} W m-2 '
              f'(table {p0_table:.2f}), flat melt at dT -4..+6: '
              + ' '.join(f'{self.flat[d]:.0f}' for d in pf.DTS))

    def S(self, slope_deg, aspect_deg):
        """All-sky shortwave [NY, days, 24] on a plane without a terrain horizon."""
        sl, asp = np.radians([slope_deg]), np.radians([aspect_deg])
        sdir = self.sun.direct(self.beam, sl, asp, np.zeros((1, 72)))[0]
        svf = float(self.ec.sky_view(sl, asp, np.zeros((1, len(self.ec.AZ))))[0])
        return self.sig * sdir[None] + self.dif * svf * self.fdif[None]

    def F(self, slope_deg, aspect_deg, dT=0.0):
        return self.pf.eb_jjas(self.Th + dT, self.S(slope_deg, aspect_deg), self.p0) / self.flat[dT]


def eb_day(ex):
    pf = ex.pf
    july = slice(181 - 151, 212 - 151)
    clear = ex.sig[:, july, 0] >= 0.8
    T = ex.Th[:, july][clear].mean(0)
    psi = np.where(T >= pf.TTIP, ex.p0 + pf.C * T, ex.p0)
    hrs = np.arange(24) + 0.5
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor=SURFACE)
    for ax in axes:
        style(ax)
        ax.set_xlim(3, 21)
        ax.set_xticks(range(4, 21, 2))
        ax.set_xlabel('hour of the day (solar time)')
    a, b = axes
    a.plot(hrs, -psi, color=INK, lw=2.2, label='sun needed before melt starts (−ψ)')
    rows = []
    for name, sl, asp, p in PLANES:
        Sab = (1 - pf.ALB) * ex.S(sl, asp)[:, july][clear].mean(0)
        melt = np.maximum(Sab + psi, 0) * 3600 / ex.oe.LF
        lab = 'flat open ground' if name == 'flat' else f'30° slope facing {name}'
        a.plot(hrs, Sab, color=PART_COL[p], lw=1.6, label=lab)
        b.plot(hrs, melt, color=PART_COL[p], lw=1.6, label=f'{lab}: {melt.sum():.1f} mm')
        rows.append((name, Sab.sum(), melt.sum()))
    a.set_ylabel('energy (W m⁻²)')
    a.set_ylim(0, 260)
    a.legend(frameon=False, fontsize=8, loc='upper left')
    a.set_title('absorbed sun against the air term', fontsize=10, color=INK)
    b.set_ylabel('melt (mm water per hour)')
    b.set_ylim(bottom=0)
    b.legend(frameon=False, fontsize=8, loc='upper left', title='melt that day', title_fontsize=8)
    b.set_title('melt = the gap, where the sun is above the line', fontsize=10, color=INK)
    fig.suptitle(f'{EX_NAME}: an average clear July day (air {T.min():+.1f} to {T.max():+.1f} °C)', fontsize=10,
                 color=INK2)
    fig.tight_layout()
    save(fig, 'ch12_eb_day')
    for name, s, m in rows:
        print(f'  clear July day, {name:5s}: absorbed sun {s * 3600 / 1e6:5.1f} MJ m-2, melt {m:4.1f} mm')
    print(f'  air {T.min():+.1f} to {T.max():+.1f} C (min at {hrs[T.argmin()]:.1f} h, max at {hrs[T.argmax()]:.1f} h);'
          f' psi0 {ex.p0:.1f} W m-2')


def eb_factor(ex):
    pf = ex.pf
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.7), facecolor=SURFACE)
    for ax in axes:
        style(ax)
        ax.axhline(1, color=INK2, lw=0.8)
    a, b, c = axes
    slopes = np.arange(0, 51, 5.0)
    print('  F at dT 0 by slope (deg):', ' '.join(f'{s:.0f}' for s in slopes))
    for name, _, asp, p in PLANES[1:]:
        f = [ex.F(s, asp) for s in slopes]
        a.plot(slopes, f, 'o-', color=PART_COL[p], ms=3.5, lw=1.6, label=f'facing {name}')
        print(f'    {name:5s}', ' '.join(f'{v:.2f}' for v in f))
    a.set_xlabel('slope (°)')
    a.set_ylabel('orientation factor F')
    a.set_title('steepness (normal warmth)', fontsize=10, color=INK)
    a.legend(frameon=False, fontsize=8, loc='lower left')
    dts = np.array(pf.DTS)
    print('  F at 30 deg by dT:', ' '.join(f'{d:+.0f}' for d in dts))
    for name, sl, asp, p in PLANES:
        f = [ex.F(sl, asp, d) for d in dts]
        b.plot(dts, f, 'o-', color=PART_COL[p], ms=3.5, lw=1.6,
               label='flat open ground' if name == 'flat' else f'30°, facing {name}')
        print(f'    {name:5s}', ' '.join(f'{v:.2f}' for v in f))
    b.set_xlabel('warmth dT of the month (°C)')
    b.set_title('warmth (30° slopes)', fontsize=10, color=INK)
    b.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    cols = [f'F{d:+.0f}' for d in pf.DTS]
    for f, ls, mk, parts, tag in (('orient_parts_v2', '-', 'o', 'GNS', 'three parts'),
                                  ('orient_parts5', '--', 's', 'GNESW', 'five parts')):
        R = pd.read_csv(W3 / f / 'parts.csv', dtype={'rgi': str})
        d = R[(R.rgi == EX_GID) & (R.band == EX_BAND)].set_index('part')
        for p in parts:
            r = d.loc[p]
            lab = f'{tag}: {p} ({r.slope:.0f}°, {r.share * 100:.0f} %)'
            c.plot(dts, r[cols].to_numpy(float), ls=ls, marker=mk, color=PART_COL[p], ms=3.5, lw=1.4, label=lab)
    c.set_xlabel('warmth dT of the month (°C)')
    c.set_title('the band\'s real parts (mean slope, share)', fontsize=10, color=INK)
    c.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    for ax in (b, c):
        ax.set_ylim(0.2, 1.25)
    fig.suptitle(EX_NAME, fontsize=10, color=INK2)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    save(fig, 'ch12_eb_factor')


def pixels():
    files = sorted((W3 / 'orient_pixels/pixels').glob('pixels_*.csv.gz'))
    P = pd.concat([pd.read_csv(f, dtype={'rgi': str}, usecols=['rgi', 'band', 'slope', 'aspect', 'F+0'])
                   for f in files], ignore_index=True)
    P = P[P.band.between(3300, 3700)]
    fig, ax = plt.subplots(figsize=(8.4, 4.6), facecolor=SURFACE)
    style(ax)
    classes = ((10, 20, LIGHT, '10–20°'), (20, 35, BLUE, '20–35°'), (35, 90, '#123a6b', 'steeper than 35°'))
    g = P[P.slope < 10]
    ax.scatter(g.aspect, g['F+0'], s=1, color=NEUTRAL, alpha=0.25, rasterized=True)
    bins = np.arange(0, 361, 15)
    mid = bins[:-1] + 7.5
    ax.plot(mid, g.groupby(pd.cut(g.aspect, bins), observed=False)['F+0'].median().to_numpy(), color=NEUTRAL, lw=2,
            label=f'under 10° (gentle): median')
    for lo, hi, col, lab in classes:
        d = P[(P.slope >= lo) & (P.slope < hi)]
        ax.scatter(d.aspect, d['F+0'], s=1, color=col, alpha=0.2, rasterized=True)
        ax.plot(mid, d.groupby(pd.cut(d.aspect, bins), observed=False)['F+0'].median().to_numpy(), color=col, lw=2,
                label=f'{lab}: median')
    for x in (90, 270):
        ax.axvline(x, color=INK, lw=1.2)
    for x in (45, 135, 225, 315):
        ax.axvline(x, color=INK2, lw=0.9, ls='--')
    ax.axhline(1, color=INK2, lw=0.8)
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 1.3)
    ax.set_xticks(range(0, 361, 45))
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N'])
    ax.set_xlabel('direction the slope faces (solid lines: three parts; dashed: five parts)')
    ax.set_ylabel('orientation factor F of the pixel')
    ax.legend(frameon=False, fontsize=8, loc='lower center', markerscale=6)
    fig.tight_layout()
    save(fig, 'ch12_pixels')
    print(f'  pixels 3300-3700 m: {len(P)} on {P.rgi.nunique()} glaciers; median F by direction and slope:')
    P['dir'] = np.array(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])[(np.mod(P.aspect + 22.5, 360) // 45).astype(int)]
    P['cls'] = pd.cut(P.slope, [0, 10, 20, 35, 90], right=False, labels=['<10', '10-20', '20-35', '>35'])
    t = P.groupby(['cls', 'dir'], observed=True)['F+0'].median().unstack()[['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']]
    print(t.round(2).to_string())
    print('  share of pixels by slope class:', (P.cls.value_counts(normalize=True).sort_index() * 100).round(0).to_dict())
    for lab, sec in (('north of three parts (270-90)', (P.aspect >= 270) | (P.aspect < 90)),
                     ('north of five parts (315-45)', (P.aspect >= 315) | (P.aspect < 45))):
        q = P[sec & (P.slope >= 10)]
        print(f'  {lab}: {len(q)} pixels, median slope {q.slope.median():.0f} deg, steeper than 35 deg '
              f'{(q.slope >= 35).mean():.2f}, median F {q["F+0"].median():.2f}')


if __name__ == '__main__':
    burial()
    air_bias()
    parts()
    if (W3 / 'orient_parts_score_v2/parts_score.json').exists():
        shade_slope()
        results()
    ex = Example()
    eb_day(ex)
    eb_factor(ex)
    if (W3 / 'orient_pixels/pixels').exists():
        pixels()
