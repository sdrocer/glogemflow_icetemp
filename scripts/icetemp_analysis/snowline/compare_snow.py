"""Snow check of a c_prec recalibration: BEFORE vs AFTER monthly forward runs, snow line against Sommer's
satellite snow lines, winter balance and ELA against GLAMOS, and the adoption rule.

Step B (plan of 2026-10-16, afternoon): the defaults (fwdB_none_fixed vs fwdB_boot7, all May-September readings,
hydrological years 1999-2025, condition (a) on the month-end values).
The monthly snow-line test (plan of 2026-10-16, evening): --rising (readings up to each glacier-year's highest
one), --years 2000 2025, --primary int (interpolated between month ends).

Usage: python compare_snow.py [OUT_DIR] [--before RUN] [--after RUN] [--changed ID ...] [--rising]
                              [--years Y0 Y1] [--primary mon|int]
"""
import argparse
import calendar
import glob
import os
import sys

import numpy as np
import pandas as pd

SNOWCAL = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal'
ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
ap.add_argument('out', nargs='?', default=os.getcwd())
ap.add_argument('--before', default=f'{SNOWCAL}/fwdB_none_fixed')
ap.add_argument('--after', default=f'{SNOWCAL}/fwdB_boot7')
ap.add_argument('--changed', nargs='+', default=['00773', '01205', '01225', '02216', '02596', '02611', '03116'])
ap.add_argument('--rising', action='store_true', help='snow-line readings up to each glacier-year\'s highest one')
ap.add_argument('--years', type=int, nargs=2, default=(1999, 2025), help='hydrological years of the readings')
ap.add_argument('--primary', choices=('mon', 'int'), default='mon', help='model value for condition (a)')
A = ap.parse_args()
OUT = A.out
RUNS = {'before': A.before, 'after': A.after}
SUB = 'monthly/CentralEurope/PAST/PAST_original'
DATA = '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/snowlines'
FIX = '/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier'
MH = '/scratch_net/vierzack06_fourth/mhuss/MassBalance_all'
Y0, NY = 1999, 27
CHANGED = A.changed
TSLA = ['00773', '00774', '01112', '01127', '01205', '01225', '01230', '02216', '02596', '02611', '02821', '03116']
MEAS_W = (1, 2, 4, 5)  # GLAMOS methods with a measured winter balance

tg = pd.read_csv(f'{DATA}/targets.csv')
NAME = {r.rgi7[-5:]: r['name'] for _, r in tg.iterrows()}
MATCH = pd.read_csv(f'{DATA}/reference/glamos_vs_geodetic.csv', dtype={'rgi7': str})


def table(run, var):
    fn = f'{RUNS[run]}/{SUB}/centraleurope_{var}_r1_CentralEurope_fair.dat'
    lines = open(fn).readlines()
    out = {}
    for ln in lines[1:]:
        p = ln.split()
        if p:
            out[p[0]] = np.array(p[1:], float)
    return out


def calib(run):
    fn = glob.glob(f'{RUNS[run]}/monthly/CentralEurope/calibration/*.dat')[0]
    return pd.read_csv(fn, sep=r'\s+', dtype={'ID': str}).set_index('ID')


def model_month(arr, year, month):
    Y = year + 1 if month >= 10 else year
    i = Y - Y0
    return arr[i, (month - 10) % 12] if 0 <= i < NY else np.nan


def model_interp(arr, t):
    nd = calendar.monthrange(t.year, t.month)[1]
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = model_month(arr, prev.year, prev.month), model_month(arr, t.year, t.month)
    return a + (b - a) * t.day / nd


def wstats(d, mcol):
    ok = np.isfinite(d[mcol]) & np.isfinite(d.sla_raw_m)
    w, e = 1 / d.sigma_m[ok], (d[mcol] - d.sla_raw_m)[ok]
    if ok.sum() == 0:
        return np.nan, np.nan, 0
    return np.sqrt(np.sum(w * e * e) / w.sum()), np.sum(w * e) / w.sum(), int(ok.sum())


# ---------------- 1. snow line ----------------
sl = {r: {k: v.reshape(NY, 12) for k, v in table(r, 'Snowline_mon').items()} for r in RUNS}
obs = []
for gid in TSLA:
    d = pd.read_csv(f'{DATA}/tsla_csv/tsla_RGI70-11.{gid}.csv', parse_dates=['date'])
    d = d[(d.qa >= 0.666) & d.month.between(5, 9) & d.hydro_year.between(*A.years)].copy()
    if A.rising:
        top = d.loc[d.groupby('hydro_year').sla_raw_m.idxmax()].set_index('hydro_year').date
        d = d[d.date <= d.hydro_year.map(top)].copy()
    d['gid'] = gid
    for r in RUNS:
        d[f'{r}_mon'] = [model_month(sl[r][gid], t.year, t.month) for t in d.date]
        d[f'{r}_int'] = [model_interp(sl[r][gid], t) for t in d.date]
    obs.append(d)
obs = pd.concat(obs, ignore_index=True)
GC, GU = f'changed{len(CHANGED)}', f'unchanged{len(TSLA) - len(CHANGED)}'
obs['group'] = np.where(obs.gid.isin(CHANGED), GC, GU)
obs.to_csv(f'{OUT}/snowline_matched_obs.csv', index=False, float_format='%.1f')


def stats_block(d):
    rec = {}
    for var in ('mon', 'int'):
        for r in RUNS:
            rm, b, n = wstats(d, f'{r}_{var}')
            rec[f'{r}_wrmse_{var}'] = rm
            rec[f'{r}_wbias_{var}'] = b
        rec[f'n'] = n
        rec[f'd_wrmse_{var}'] = rec[f'after_wrmse_{var}'] - rec[f'before_wrmse_{var}']
    return rec


rows = []
for gid, d in obs.groupby('gid'):
    rows.append(dict(rgi7=gid, name=NAME[gid], changed=gid in CHANGED,
                     yrs=f'{d.hydro_year.min()}-{d.hydro_year.max()}', **stats_block(d)))
for grp, d in obs.groupby('group'):
    rows.append(dict(rgi7=f'POOLED_{grp}', name='pooled (obs-weighted 1/sigma)', changed=grp == GC,
                     yrs=f'{d.hydro_year.min()}-{d.hydro_year.max()}', **stats_block(d)))
SLG = pd.DataFrame(rows)
# glacier-mean of per-glacier scores (each glacier weighted equally) as a sensitivity
for grp, ch in ((GC, True), (GU, False)):
    s = SLG[(SLG.changed == ch) & ~SLG.rgi7.str.startswith('POOLED')]
    rec = dict(rgi7=f'GLACIERMEAN_{grp}', name='mean of per-glacier scores', changed=ch, yrs='')
    for c in s.columns:
        if c.startswith(('before_', 'after_', 'd_')):
            rec[c] = s[c].mean()
    rec['n'] = s.n.sum()
    rows.append(rec)
SLG = pd.DataFrame(rows)
SLG.to_csv(f'{OUT}/snowline_per_glacier.csv', index=False, float_format='%.1f')

rows = []
for grp, dg in obs.groupby('group'):
    for m, d in dg.groupby('month'):
        rows.append(dict(group=grp, month=m, **stats_block(d)))
SLM = pd.DataFrame(rows)
SLM.to_csv(f'{OUT}/snowline_by_month.csv', index=False, float_format='%.1f')

# per glacier and month (changed only)
rows = []
for (gid, m), d in obs[obs.gid.isin(CHANGED)].groupby(['gid', 'month']):
    rows.append(dict(rgi7=gid, name=NAME[gid], month=m, **stats_block(d)))
pd.DataFrame(rows).to_csv(f'{OUT}/snowline_by_glacier_month_{GC}.csv', index=False, float_format='%.1f')

# sensitivity: drop inferred full-snow / snow-free observations
rows = []
for grp, d in obs[obs.flag.isna()].groupby('group'):
    rows.append(dict(group=grp, variant='no_inferred_flags', **stats_block(d)))
pd.DataFrame(rows).to_csv(f'{OUT}/snowline_pooled_noflag.csv', index=False, float_format='%.1f')


# ---------------- 2/3. GLAMOS winter balance, ELA ----------------
def read_fix(name):
    lines = open(f'{FIX}/{name}_fix.dat', encoding='latin-1').readlines()
    h = [s.strip() for s in lines[0].split(';')]
    n, zlo, zhi = int(h[4]), float(h[5]), float(h[6])
    edges = np.linspace(zlo, zhi, n + 1)
    recs = {}
    for l in lines:
        if l.startswith('#') or not l.strip():
            continue
        p = l.split()
        f = lambda s: np.nan if s.lower() == 'nan' else float(s.rstrip('.'))
        b = np.array([f(x) for x in p[12:]])
        bw, nw, ba, na, A = (b[i * n:(i + 1) * n] for i in range(5))
        recs[int(p[4][:4])] = dict(meth=int(p[0]), Bw=f(p[5]) / 1000, Ba=f(p[6]) / 1000, ELA=f(p[7]),
                                   A=f(p[9]), bw=bw / 1000, nw=nw, A_eb=A)
    return edges, recs


def read_mh(names):
    """mhuss bands_<yr>.dat summed over several glaciers: lo, hi, area, fxd_ann, meas_ann, fxd_win, meas_win."""
    sets = []
    for nm in names:
        d = f'{MH}/{nm}/massbalance/results'
        s = {}
        for fn in sorted(glob.glob(f'{d}/bands_*.dat')):
            rows = [[float(x) for x in l.split()[:7]] for l in open(fn).readlines()[1:] if len(l.split()) >= 7]
            s[int(os.path.basename(fn)[6:10])] = np.unique(np.array(rows), axis=0)
        sets.append(s)
    lo = min(a[:, 0].min() for s in sets for a in s.values())
    hi = max(a[:, 1].max() for s in sets for a in s.values())
    a0 = next(iter(sets[0].values()))
    w = a0[0, 1] - a0[0, 0]
    edges = np.arange(lo, hi + w / 2, w)
    n = len(edges) - 1
    recs = {}
    for yr in sorted(set.intersection(*[set(s) for s in sets])):
        A = np.zeros(n); W = np.zeros(n)
        for s in sets:
            for r in s[yr]:
                i = int(round((r[0] - lo) / w))
                if 0 <= i < n and r[2] > 0:
                    A[i] += r[2]; W[i] += r[2] * r[5]
        with np.errstate(all='ignore'):
            bw = np.where(A > 0, W / A, np.nan) / 1000
        recs[yr] = dict(meth=-1, Bw=np.nansum(bw * A) / A.sum(), ELA=np.nan, A=A.sum(), bw=bw, nw=np.zeros(n), A_eb=A)
    return edges, recs


def mband(run, rid):
    rows = [l.split() for l in open(f'{RUNS[run]}/{SUB}/mb_elevation/bwelev_{rid}.dat').readlines()[1:] if l.strip()]
    z = np.array([float(r[0]) for r in rows])
    v = np.array([[float(x) for x in r[1:]] for r in rows])
    v[v <= -98] = np.nan
    return z, v


def on_bins(z, v, edges):
    out = np.full((len(edges) - 1, v.shape[1]), np.nan)
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1])
        if m.any():
            with np.errstate(all='ignore'):
                out[i] = np.nanmean(v[m], axis=0)
    return out


WB = {r: table(r, 'Winter_balance_sfc') for r in RUNS}
EL = {r: table(r, 'ELA') for r in RUNS}
AR = {r: table(r, 'Area') for r in RUNS}

# rgi7, label, source, outline ratio, flag, counts for condition (b)
ratio = {r.short: r.area_ratio for _, r in MATCH.iterrows()}
SRC = [('02596', 'fix:aletsch', True), ('02216', 'fix:pers', True),
       ('01112', 'fix:hohsaas', False), ('01309', 'fix:sexrouge', False), ('01039', 'fix:tortin', False),
       ('01150', 'fix:alphubel', False), ('02103', 'fix:corvatsch', False),
       ('02216', 'mh:morteratsch+pers', False), ('02611', 'mh:fiescher', False), ('01225', 'mh:gorner', False)]

bw_rows, ela_rows, band_rows, band_prof = [], [], [], []
for rid, src, cond_b in SRC:
    kind, nm = src.split(':')
    if kind == 'fix':
        edges, recs = read_fix(nm)
        rat = ratio.get(nm, np.nan)
    else:
        edges, recs = read_mh(nm.split('+'))
        rat = np.nan
    yrs_all = [y for y in recs if 2000 <= y <= 2025 and np.isfinite(recs[y]['Bw'])]
    if kind == 'fix':
        yrs_meas = [y for y in yrs_all if recs[y]['meth'] in MEAS_W]
    elif nm == 'morteratsch+pers':
        yrs_meas = [y for y in yrs_all if 2021 <= y <= 2025]  # winter stakes on both (peg_*_w.dat)
    else:
        yrs_meas = []  # no winter stakes in 2000-2025 (fiescher: first winter survey 2026; gorner: none)
    A_src = np.nanmean([recs[y]['A'] for y in yrs_all]) if yrs_all else np.nan
    if kind == 'mh':
        rat = A_src / AR['before'][rid][2003 - Y0 + 0] if rid in AR['before'] else np.nan
    flag = '' if (np.isfinite(rat) and 0.85 <= rat <= 1.15) else f'outline ratio {rat:.2f}'
    if kind == 'mh':
        flag = ('mhuss evaluation, not GLAMOS-published; ' + ('winter stakes 2021-2025' if yrs_meas else
                'NO winter stakes 2000-2025 -> winter is model-inferred, info only') + '; ' + flag)
    for sel, yrs in (('measured', yrs_meas), ('all', yrs_all)):
        if not yrs:
            continue
        o = np.mean([recs[y]['Bw'] for y in yrs])
        rec = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(yrs),
                   yrs=f'{min(yrs)}-{max(yrs)}', area_obs_km2=A_src, outline_ratio=rat, flag=flag.strip('; '),
                   Bw_obs=o)
        for r in RUNS:
            m = np.mean([WB[r][rid][y - Y0] for y in yrs])
            rec[f'Bw_{r}'] = m
            rec[f'absdiff_{r}'] = abs(m - o)
        rec['cond_b_pass'] = rec['absdiff_after'] <= rec['absdiff_before'] + 1e-9
        rec['counts_for_cond_b'] = cond_b and sel == 'measured' and rid in CHANGED
        bw_rows.append(rec)
        # ELA
        ey = [y for y in yrs if np.isfinite(recs[y]['ELA'])]
        if ey:
            o = np.mean([recs[y]['ELA'] for y in ey])
            rec2 = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(ey),
                        yrs=f'{min(ey)}-{max(ey)}', flag=flag.strip('; '), ELA_obs=o)
            for r in RUNS:
                rec2[f'ELA_{r}'] = np.mean([EL[r][rid][y - Y0] for y in ey])
                rec2[f'bias_{r}'] = rec2[f'ELA_{r}'] - o
            ela_rows.append(rec2)
        # bands
        ob = np.array([recs[y]['bw'] for y in yrs]).T
        Aw = np.array([recs[y]['A_eb'] for y in yrs]).T
        nwt = np.array([recs[y]['nw'] for y in yrs]).T
        zc = 0.5 * (edges[1:] + edges[:-1])
        rec3 = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(yrs),
                    yrs=f'{min(yrs)}-{max(yrs)}', bin_m=edges[1] - edges[0], flag=flag.strip('; '))
        for r in RUNS:
            z, v = mband(r, rid)
            mb = on_bins(z, v, edges)[:, [y - Y0 for y in yrs]]
            ok = np.isfinite(ob) & np.isfinite(mb) & (Aw > 0)
            w = np.where(ok, Aw, 0)
            e = np.where(ok, mb - ob, 0)
            rec3[f'n_bandyears'] = int(ok.sum())
            rec3[f'bias_{r}'] = np.sum(w * e) / w.sum()
            rec3[f'rmse_{r}'] = np.sqrt(np.sum(w * e * e) / w.sum())
            okp = ok & (nwt > 0)
            if okp.sum():
                wp = np.where(okp, Aw, 0)
                rec3['n_bandyears_with_pts'] = int(okp.sum())
                rec3[f'bias_ptbins_{r}'] = np.sum(wp * e) / wp.sum()
                rec3[f'rmse_ptbins_{r}'] = np.sqrt(np.sum(wp * e * e) / wp.sum())
            with np.errstate(all='ignore'):
                mprof = np.nanmean(np.where(ok, mb, np.nan), axis=1)
                oprof = np.nanmean(np.where(ok, ob, np.nan), axis=1)
            for i in np.where(np.isfinite(oprof))[0]:
                band_prof.append(dict(rgi7=rid, source=src, years_sel=sel, run=r, z_lo=edges[i], z_hi=edges[i + 1],
                                      A_obs=np.nanmean(Aw[i]), n_winter_pts=int(nwt[i].sum()), n_years=int(ok[i].sum()),
                                      bw_obs=oprof[i], bw_model=mprof[i]))
        band_rows.append(rec3)

BW = pd.DataFrame(bw_rows); BW.to_csv(f'{OUT}/glamos_Bw_glacierwide.csv', index=False, float_format='%.3f')
ELA = pd.DataFrame(ela_rows); ELA.to_csv(f'{OUT}/glamos_ELA.csv', index=False, float_format='%.0f')
BB = pd.DataFrame(band_rows); BB.to_csv(f'{OUT}/glamos_Bw_bands_summary.csv', index=False, float_format='%.3f')
BP = pd.DataFrame(band_prof)
BP = BP.pivot_table(index=['rgi7', 'source', 'years_sel', 'z_lo', 'z_hi', 'A_obs', 'n_winter_pts', 'bw_obs'],
                    columns='run', values='bw_model').reset_index()
BP.to_csv(f'{OUT}/glamos_Bw_band_profiles.csv', index=False, float_format='%.3f')

# ---------------- 4. sanity: unchanged glaciers ----------------
SAN_N = len([g for g in table('before', 'Area') if g not in CHANGED])
rows = []
VARS = ['Snowline_mon', 'Winter_balance_sfc', 'Annual_Balance_sfc', 'ELA', 'Accumulation_sfc', 'Snowmelt_sfc',
        'Icemelt_sfc', 'Volume', 'Area', 'AAR', 'Precipitation_mon', 'Melt_mon']
for var in VARS:
    a, b = table('before', var), table('after', var)
    for gid in a:
        if gid in CHANGED:
            continue
        d = b[gid] - a[gid]
        i = int(np.nanargmax(np.abs(d))) if np.isfinite(d).any() else 0
        rows.append(dict(var=var, rgi7=gid, name=NAME.get(gid, ''), max_absdiff=np.nanmax(np.abs(d)),
                         mean_absdiff=np.nanmean(np.abs(d)), mean_before=np.nanmean(a[gid]), idx_max=i))
ca, cb = calib('before'), calib('after')
for gid in ca.index:
    if gid in CHANGED:
        continue
    for c in ('DDFsnow', 'Cprec', 'T_off', 'Ba', 'Bw'):
        rows.append(dict(var=f'calib_{c}', rgi7=gid, name=NAME.get(gid, ''), max_absdiff=abs(cb.loc[gid, c] - ca.loc[gid, c]),
                         mean_absdiff=abs(cb.loc[gid, c] - ca.loc[gid, c]), mean_before=ca.loc[gid, c], idx_max=0))
SAN = pd.DataFrame(rows); SAN.to_csv(f'{OUT}/sanity_unchanged{SAN_N}.csv', index=False, float_format='%.4g')
cal = ca[['DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'Ba']].join(cb[['DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'Ba']],
                                                                    lsuffix='_before', rsuffix='_after')
cal.insert(0, 'name', [NAME.get(g, '') for g in cal.index])
cal.insert(1, 'changed', [g in CHANGED for g in cal.index])
cal.to_csv(f'{OUT}/calibration_before_after.csv', float_format='%.3f')

# ---------------- 5. adoption rule ----------------
pc = SLG[SLG.rgi7 == f'POOLED_{GC}'].iloc[0]
cond_a = pc[f'after_wrmse_{A.primary}'] < pc[f'before_wrmse_{A.primary}']
# year-block bootstrap of the pooled wRMSE difference (resample hydrological years, all 7 glaciers together)
rng = np.random.default_rng(0)
d7 = obs[obs.group == GC]
yrs = np.sort(d7.hydro_year.unique())
byy = {y: d7[d7.hydro_year == y] for y in yrs}
boot = []
for _ in range(2000):
    dd = pd.concat([byy[y] for y in rng.choice(yrs, len(yrs))])
    boot.append([wstats(dd, 'after_mon')[0] - wstats(dd, 'before_mon')[0], wstats(dd, 'after_int')[0] - wstats(dd, 'before_int')[0]])
boot = np.array(boot)
BOOT = pd.DataFrame(dict(variant=['model_month', 'model_interp'],
                         d_wrmse=[pc.d_wrmse_mon, pc.d_wrmse_int],
                         p05=np.percentile(boot, 5, axis=0), p50=np.percentile(boot, 50, axis=0),
                         p95=np.percentile(boot, 95, axis=0), frac_improved=(boot < 0).mean(axis=0)))
BOOT.to_csv(f'{OUT}/snowline_pooled_bootstrap_{GC}.csv', index=False, float_format='%.2f')
print(BOOT.round(2))
cb_rows = BW[BW.counts_for_cond_b]
cond_b = bool(cb_rows.cond_b_pass.all())
with open(f'{OUT}/verdict.txt', 'w') as f:
    if A.primary == 'mon' and not A.rising:
        f.write(f'(a) pooled wRMSE {GC} May-Sep (model_month): before {pc.before_wrmse_mon:.1f} m, after {pc.after_wrmse_mon:.1f} m -> {cond_a}\n')
        f.write(f'    sensitivity model_interp: before {pc.before_wrmse_int:.1f}, after {pc.after_wrmse_int:.1f}\n')
    else:
        sel = 'up to each year\'s highest reading' if A.rising else 'May-Sep'
        f.write(f'(a) pooled wRMSE {GC}, {sel}, {A.years[0]}-{A.years[1]} (model_{"interp" if A.primary == "int" else "month"}): '
                f'before {pc[f"before_wrmse_{A.primary}"]:.1f} m, after {pc[f"after_wrmse_{A.primary}"]:.1f} m -> {cond_a}\n')
        other = 'mon' if A.primary == 'int' else 'int'
        f.write(f'    sensitivity model_{"month" if other == "mon" else "interp"}: before {pc[f"before_wrmse_{other}"]:.1f}, '
                f'after {pc[f"after_wrmse_{other}"]:.1f}\n')
    for _, r in cb_rows.iterrows():
        f.write(f'(b) {r.rgi7} {r["name"]} [{r.source}, {r.yrs}, n={r.n_years}] obs {r.Bw_obs:.3f} before {r.Bw_before:.3f} '
                f'after {r.Bw_after:.3f}: |d| {r.absdiff_before:.3f} -> {r.absdiff_after:.3f} pass={r.cond_b_pass} {r.flag}\n')
    f.write(f'(b) all -> {cond_b}\n')

pd.set_option('display.width', 250, 'display.max_columns', 40)
print(SLG.round(1).to_string(index=False))
print(SLM.round(1).to_string(index=False))
print(BW.round(3).to_string(index=False))
print(ELA.round(0).to_string(index=False))
print(BB.round(3).to_string(index=False))
print(SAN.sort_values('max_absdiff', ascending=False).groupby('var').head(2).round(4).to_string(index=False))
print(cal.round(3).to_string())
print(open(f'{OUT}/verdict.txt').read())
