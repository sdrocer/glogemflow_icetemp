"""Week-1 diagnosis: the 2950-3550 m zone and the frozen bed, from existing runs only.

A. Zone entities (model band 2950-3550 m): measured and modelled (today, F, H, FH) mean below
   15 m, three-class labels, miss type, and per-entity diagnostics for three candidate causes:
   (1) winter snow too thick, (2) meltwater too deep, (3) heat kept from the firn years.
   Plus every profile of the zone measured against the model at its own dates (Tete Rousse
   repeats 2010 / 2015 / 2022-23 among them).
B. Alpine holes that reached the bed: measured deepest reading against the model's bottom node
   at the hole's band, frozen or temperate by the same rule, per run.
C. Whether temp_bedrock_<id>.dat can give a glacier-level bed type, with a check of what it holds.

Conventions (checked in the IDL source):
- temp_ID rows: year = hydrological year label (Oct of year-1 .. Sep of year), month = calendar
  month, profile at the end of the month (runner.match_model_to_obs).
- firnice_snowmay: column k = snow on 1 May of calendar year 1972+k [m w.e.]; snow left from the
  previous balance year is removed each year (glogem.pro:524), so it is Oct-Apr accumulation
  minus melt, comparable with a GLAMOS fixed-date winter balance (1 Oct - 30 Apr).
- firnice_firnline / firnice_wetdepth: year = hydrological year label.
- A band is firn in year Y when band >= firnline_th[Y] (lowest thermal firn band). The mask is not
  logged per band, so isolated ice gaps above the firn line are missed (n_firn_to_ice counts them).
- On firn bands the seasonal snow is thermally inert (firnice_temperature_model.pro: snow and
  firn share fit_dens); it insulates only on ice bands.
- The column at a borehole band runs to max(band thickness, first node at or below the deepest
  observation in that band) (firnice_temperature_model.pro:95-108, since 2026-09-24). The bottom
  valid node of temp_ID is therefore the model's bed only where the band is thicker than the
  deepest observation.

Usage: PYTHONPATH=src python scripts/icetemp_analysis/zone3000_bed_check.py [out_dir]
"""
import json
import pickle
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from forcing_check import INV_HF, model_band, read_inventory, read_lookup  # noqa: E402

from icetemp.calibration.runner import match_model_to_obs, model_at_depth  # noqa: E402
from icetemp.calibration.thermal_structure import (  # noqa: E402
    DEFAULT_TOLERANCE, borehole_tolerance, classify, classify_point, load_borehole_table,
    pressure_melting_point)
from icetemp.io import read_profile_file, read_summary_file  # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1
           else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week1/zone3000_bed')
POINTS = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/thermal_forcing/score/points.csv')
CACHE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/cache/datahandler.pkl')
RUNS = {'today': 'ce25_pred_w0', 'F': 'ce25_F_w0', 'H': 'ce25_H_w0', 'FH': 'ce25_FH_w0'}
SUB = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature'
BANDS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope')
BANDS_HF = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands_HF/centraleurope')
HF_AREA = 0.25               # km2: read_hypsometryfile.pro reads bands_HF below this (RGI7)
GLENGLAT = REPO / 'glenglat/data'
GLIMS_LOOKUP = REPO / 'data/rgi7_centraleurope_glims_lookup.csv'
GLAMOS = Path('/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier')

# ---- selections and thresholds, fixed before any result was computed -------------------------
ZONE = (2950, 3550)          # model band [m a.s.l.], inclusive
DEEP = 15.0                  # points strictly deeper than this [m] (forcing_runs_score DEEP)
NEAR_PMP = 1.0               # miss with |mean error below 15 m| < this [C] = near-PMP flip
SNOW_YEARS = 5               # 1 May snow: the 5 Mays before the reading
WET_YEARS = 10               # wet depth and firn status: the 10 hydrological years before it
WET_RUN = 'H'                # only H/FH/P log the wet depth
MELT_DONE_MONTH = 8          # a reading in Aug/Sep counts its own hydrological year as melted
GLAMOS_MAP = {'02103': 'corvatsch', '01112': 'hohsaas', '02596': 'aletsch'}
SNOW_SUPPORT, SNOW_REFUTE = 1.2, 1.0   # model / GLAMOS band winter balance
NEAR_ZERO = 0.3              # cause 3: modelled deep mean at the reading within this of PMP [C]
MEMORY_YEARS = 30            # cause 3: firn lost at most this many years before the reading
ALPS = dict(lat=(43.5, 48.5), lon=(5.0, 17.0))
EQ_KEEP = {'true', 'estimated', ''}    # equilibrium states used (as DataHandler); 'false' dropped
WINDOW = (197110, 202509)    # run months (hydrological years 1972-2025)
BED_GAP_FLAG = 5.0           # flag holes whose deepest reading is > this above the hole depth [m]
LOGGER_PROFILES = 20         # holes with more profiles are loggers: pooled per calendar month
BEDROCK_CAP = -0.347         # temp_bedrock floor value (PMP of the 519 m slot; see part C)
TYPE_YEARS = (1980, 2000, 2010, 2025)
NODES = np.array([0.] + list(range(1, 10)) + list(range(14, 60, 5)) + list(range(79, 520, 20)))


def run_dir(k):
    return REPO / 'data/bayescal' / RUNS[k] / SUB


def ym_add(ym, k):
    i = (ym // 100) * 12 + ym % 100 - 1 + k
    return (i // 12) * 100 + i % 12 + 1


def hyd_year(ym):
    return ym // 100 + (ym % 100 >= 10)


def reading_months(ym0, ym1, frac):
    """data._reading_months: one month at its fraction, or a range at mid-month, inside the run."""
    if ym0 == ym1:
        return np.array([ym0]), np.array([1.0]), np.array([frac])
    m = month_range(max(ym0, WINDOW[0]), min(ym1, WINDOW[1]))
    if not m:
        return np.array([ym0]), np.array([1.0]), np.array([0.5])
    return np.array(m), np.full(len(m), 1.0 / len(m)), np.full(len(m), 0.5)


def month_range(a, b):
    out = []
    while a <= b:
        out.append(a)
        a = ym_add(a, 1)
    return out


def date_months(d0, d1):
    """DataHandler's dating of a profile: midpoint month for up to 31 days, else the month range."""
    t0 = pd.to_datetime(d0, errors='coerce')
    t1 = pd.to_datetime(d1, errors='coerce')
    if pd.isna(t1):
        t1 = t0
    if pd.isna(t0):
        return None
    if (t1 - t0).days > 31:
        return t0.year * 100 + t0.month, t1.year * 100 + t1.month, 0.5
    tm = t0 + (t1 - t0) / 2
    return tm.year * 100 + tm.month, tm.year * 100 + tm.month, tm.day / tm.days_in_month


# ---- readers ---------------------------------------------------------------------------------
_inv = {}


def band_table(gid):
    """The hypsometry the run read: bands_HF for glaciers under HF_AREA km2 (files_HF area)."""
    if not _inv:
        ids, v = read_inventory(INV_HF)
        _inv.update(zip(ids, v[:, 2]))
    src = BANDS_HF if _inv.get(gid, 1e9) < HF_AREA and (BANDS_HF / f'{gid}.dat').exists() else BANDS
    rows = []
    for ln in (src / f'{gid}.dat').read_text().splitlines()[5:]:
        a = ln.split()
        if len(a) >= 5:
            rows.append(dict(band=int(a[1]) + 5, area=float(a[3]), thick=float(a[4])))
    return pd.DataFrame(rows).set_index('band')


def read_firnline(k, gid):
    a = np.loadtxt(run_dir(k) / f'firnice_firnline_{gid}.dat', comments='#')
    return pd.DataFrame(a, columns=['year', 'ela', 'aar', 'firnline', 'n_f2i', 'n_i2f']).astype(
        {'year': int}).set_index('year')


def read_snowmay(gid):
    rows = {}
    for ln in (run_dir('today') / f'firnice_snowmay_{gid}.dat').read_text().splitlines():
        if ln.strip() and not ln.startswith('#'):
            a = ln.split()
            rows[int(a[0])] = np.array(a[1:], dtype=float)
    return pd.DataFrame(rows, index=range(1972, 1972 + len(next(iter(rows.values()))))).T


def read_wetdepth(k, gid):
    p = run_dir(k) / f'firnice_wetdepth_{gid}.dat'
    if not p.exists():
        return None
    a = np.loadtxt(p, comments='#', ndmin=2)
    return pd.DataFrame(a, columns=['elev', 'year', 'wet', 'lost']).astype({'elev': int, 'year': int})


_prof_cache = {}


def profile(k, gid, pid):
    key = (k, gid, pid)
    if key not in _prof_cache:
        p = run_dir(k) / f'{pid}_{gid}.dat'
        elev, depths, df = read_profile_file(p)
        _prof_cache[key] = (elev, np.asarray(depths, float), df)
    return _prof_cache[key]


def glamos_bands(name):
    p = GLAMOS / f'{name}_fix.dat'
    lines = p.read_text().splitlines()
    h = [s.strip() for s in lines[0].split(';')]
    n, lo, hi = int(h[4]), float(h[5]), float(h[6])
    rows = []
    for ln in lines:
        if ln.startswith('#') or not ln.strip():
            continue
        a = ln.split()
        rows.append(dict(year=int(a[4][:4]), date_s=a[3], bw=float(a[5]) / 1000, ela=float(a[7]),
                         hmax=float(a[11]),
                         bands=np.array(a[12:12 + n], dtype=float) / 1000,
                         area=np.array(a[12 + 4 * n:12 + 5 * n], dtype=float)))
    return dict(n=n, lo=lo, hi=hi, w=(hi - lo) / n, rows=rows, path=str(p))


def model_series(df, d):
    cols = [c for c in df.columns if c not in ('year', 'month')]
    return model_at_depth(np.asarray(cols, float), df[cols].to_numpy(float), d)


def bottom_node(df):
    cols = [c for c in df.columns if c not in ('year', 'month')]
    T = df[cols].to_numpy(float)
    res = np.flatnonzero(np.isfinite(T).any(axis=0))
    return (float(cols[res.max()]), T[:, res.max()]) if res.size else (np.nan, None)


def model_at(df, depths, months_list):
    """match_model_to_obs for given depths, each with its (months, weights, fractions)."""
    g = SimpleNamespace(depths=np.asarray(depths, float), obs_months=months_list)
    _, v = match_model_to_obs(df, g)
    return v[:len(depths)]


# ---- part A ----------------------------------------------------------------------------------
def firn_status(fl, band):
    return {int(y): bool(band >= r.firnline) for y, r in fl.iterrows()}


def firn_loss_year(status):
    ys = sorted(status)
    if all(not status[y] for y in ys):
        return 'ice throughout'
    if status[ys[-1]]:
        return 'firn at end'
    for y in ys:
        if all(not status[z] for z in ys if z >= y):
            return y
    return None


def part_a(P, gl, raw, bdf, lookup):
    by_name = {g.glacier_name: g for g in gl}
    zone = P[(P.band >= ZONE[0]) & (P.band <= ZONE[1])]
    ent_rows, prof_rows, series_rows, glamos_rows = [], [], [], []
    glamos_cache = {}
    for ent, sub in zone.groupby('entity'):
        g = by_name[ent]
        gid = g.glacier_id
        tol = borehole_tolerance(g.borehole_id, bdf)
        snow = read_snowmay(gid)
        band, src = model_band(g, lookup, snow.index.to_numpy())
        pid = src if src.startswith('temp_ID') else None
        deep = sub[sub.depth > DEEP]
        row = dict(entity=ent, rgi=gid, band_entity=int(g.elevation_band), band_model=band,
                   profile=pid, tolerance=tol, n_points=len(sub), n_deep=len(deep),
                   labels='/'.join(sorted(set(sub.label))))
        use = deep if len(deep) else sub
        row['set'] = 'below 15 m' if len(deep) else 'all depths (none below 15 m)'
        row['obs'] = float(use.obs.mean())
        row['class_obs'] = classify(use.obs, use.depth, tol)
        for k in RUNS:
            row[k] = float(use[k].mean())
            row[f'class_{k}'] = classify(use[k], use.depth, tol)
            row[f'err_{k}'] = row[k] - row['obs']
        if len(deep):
            hit = row['class_obs'] == row['class_today']
            row['miss'] = 'hit' if hit else ('near-PMP flip' if abs(row['err_today']) < NEAR_PMP
                                             else 'clear miss')
        else:
            row['miss'] = 'no point below 15 m'
        # profiles behind the entity, rebuilt from glenglat as DataHandler selects them
        rd = raw[(raw.glacier_name == g.base_glacier_name) & (raw.band == int(g.elevation_band))]
        rdeep = rd[rd.depth > DEEP] if len(deep) else rd
        row['n_readings_raw'] = len(rd)
        row['n_readings_cache'] = int(np.sum(g.n_readings))
        read_ym = sorted(set(int(v) for v in rdeep.ym1))
        row['reading_months'] = ','.join(map(str, read_ym))
        fl = read_firnline('today', gid)
        status = firn_status(fl, band)
        loss = firn_loss_year(status)
        row['firn_loss_year'] = loss
        row['firn_years'] = int(sum(status.values()))
        for dec in (1970, 1980, 1990, 2000, 2010, 2020):
            ys = [y for y in fl.index if dec <= y < dec + 10]
            row[f'firnline_{dec}s'] = float(fl.loc[ys, 'firnline'].mean())
        # (1) snow, (2) wet depth, firn status before each reading
        wd = read_wetdepth(WET_RUN, gid)
        wdb = wd[wd.elev == band].set_index('year') if wd is not None else None
        d_obs = float(use.depth.min())
        snow_pre, n_ice5, n_firn10, n_wet10, wet_max10, read_hy, glam = [], 0, 0, 0, 0.0, [], []
        for ym in read_ym:
            y, mo = ym // 100, ym % 100
            mays = list(range(y - SNOW_YEARS + 1, y + 1)) if mo >= 5 else list(range(y - SNOW_YEARS, y))
            mays = [m for m in mays if m in snow.columns]
            snow_pre.append(float(snow.loc[band, mays].mean()))
            hy = hyd_year(ym)
            read_hy.append(hy)
            last = hy if MELT_DONE_MONTH <= mo <= 9 else hy - 1
            win = [z for z in range(last - WET_YEARS + 1, last + 1) if z in status]
            n_ice5 += sum(1 for m in mays if m in status and not status[m])
            n_firn10 += sum(1 for z in win if status[z])
            if wdb is not None and win:
                w = wdb.reindex(win)
                n_wet10 += int(((w.wet >= d_obs) & np.array([status[z] for z in win])).sum())
                wet_max10 = max(wet_max10, float(np.nanmax(w.wet.to_numpy())) if w.wet.notna().any() else 0.0)
            if gid in GLAMOS_MAP:
                glam.append(mays)
        row['snow_may_pre5'] = float(np.mean(snow_pre))
        row['snow_may_mean_run'] = float(snow.loc[band].mean())
        row['n_mays_band_ice'] = n_ice5
        row['n_mays_total'] = SNOW_YEARS * len(read_ym)
        row['n_years_firn_pre10'] = n_firn10
        row['n_years_wet_to_obs_pre10'] = n_wet10
        row['wet_max_pre10'] = wet_max10
        row['d_obs_min'] = d_obs
        row['firn_at_reading'] = '/'.join('firn' if status.get(h) else 'ice' for h in read_hy)
        # the ice run that ends at the first reading: its first year is when the band lost its firn
        loss_rel = None
        if read_hy:
            h0 = min(read_hy)
            if h0 in status and not status[h0]:
                y = h0
                while y - 1 in status and not status[y - 1]:
                    y -= 1
                loss_rel = y if y - 1 in status else 'ice since 1972'
            elif h0 in status:
                loss_rel = 'firn at reading'
        row['firn_lost_before_reading'] = loss_rel
        row['wet_years_ge15_run'] = int((wdb.wet >= 15).sum()) if wdb is not None else None
        row['wet_last_year_ge_obs'] = (int(wdb[wdb.wet >= d_obs].index.max())
                                       if wdb is not None and (wdb.wet >= d_obs).any() else None)
        # GLAMOS band winter balance at the model band, overlapping years
        if gid in GLAMOS_MAP:
            name = GLAMOS_MAP[gid]
            gb = glamos_cache.setdefault(name, glamos_bands(name))
            ib = int((band - gb['lo']) // gb['w'])
            areas = band_table(gid)
            pairs = []
            for r in gb['rows']:
                if 0 <= ib < gb['n'] and np.isfinite(r['bands'][ib]) and r['year'] in snow.columns:
                    mw = float(snow.loc[band, r['year']])
                    gw_mod = float((snow[r['year']] * areas.reindex(snow.index).area).sum()
                                   / areas.reindex(snow.index).area.sum())
                    pairs.append((r['year'], mw, r['bands'][ib], gw_mod, r['bw']))
                    glamos_rows.append(dict(entity=ent, rgi=gid, glamos=name, year=r['year'],
                                            band_model=band, glamos_bin=f"{gb['lo'] + ib * gb['w']:.0f}-"
                                            f"{gb['lo'] + (ib + 1) * gb['w']:.0f}",
                                            model_snow_may=mw, glamos_bw_bin=r['bands'][ib],
                                            model_glacier_wide=gw_mod, glamos_bw_glacier=r['bw']))
            if pairs:
                a = np.array([p[1:] for p in pairs])
                pre = [p for p in pairs if any(p[0] in m for m in glam)]
                row.update(glamos=name, glamos_years=f'{pairs[0][0]}-{pairs[-1][0]} ({len(pairs)})',
                           glamos_overlap_pre5=len(pre), model_band_mean=float(a[:, 0].mean()),
                           glamos_band_mean=float(a[:, 1].mean()),
                           ratio_band=float(a[:, 0].mean() / a[:, 1].mean()),
                           ratio_glacier=float(a[:, 2].mean() / a[:, 3].mean()))
        # (3) temperature at 15/20 m per hydrological year, every run
        if pid:
            for k in RUNS:
                _, _, df = profile(k, gid, pid)
                for d in (15.0, 20.0):
                    s = pd.Series(model_series(df, d), index=df.year.to_numpy()).groupby(level=0).mean()
                    for y, v in s.items():
                        series_rows.append(dict(entity=ent, rgi=gid, band=band, run=k, depth=d, year=int(y),
                                                T=float(v), firn=status.get(int(y))))
            S = pd.DataFrame([r for r in series_rows if r['entity'] == ent and r['run'] == 'today'
                              and r['depth'] == 20.0]).set_index('year')['T']
            ref = loss_rel if isinstance(loss_rel, int) else (loss if isinstance(loss, int) else None)
            row['T20_ref_year'] = ref
            if ref is not None:
                for off in (-1, 2, 5, 10, 15):
                    row[f'T20_loss{off:+d}'] = float(S.get(ref + off, np.nan))
            row['T20_at_reading'] = float(np.nanmean([S.get(h, np.nan) for h in read_hy]))
            row['T20_1972'] = float(S.get(1972, np.nan))
            # every profile against the model at its own dates
            for (bh, unit), rp in rdeep.groupby(['borehole_id', 'unit']):
                rp = rp.groupby(['depth', 'ym0', 'ym1', 'frac'], as_index=False).agg(
                    temperature=('temperature', 'mean'), label=('label', 'first'),
                    date_min=('date_min', 'min'), date_max=('date_max', 'max'))
                months = [reading_months(int(a), int(b), float(f))
                          for a, b, f in zip(rp.ym0, rp.ym1, rp.frac)]
                rec = dict(entity=ent, rgi=gid, band=band, borehole_id=bh, unit=unit,
                           label=rp.label.iloc[0], date=rp.date_min.min(), date_end=rp.date_max.max(),
                           n=len(rp), depth_min=float(rp.depth.min()), depth_max=float(rp.depth.max()),
                           obs=float(rp.temperature.mean()),
                           class_obs=classify(rp.temperature, rp.depth, tol))
                for k in RUNS:
                    _, _, df = profile(k, gid, pid)
                    v = model_at(df, rp.depth.to_numpy(), months)
                    rec[k] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
                    rec[f'n_model_{k}'] = int(np.isfinite(v).sum())
                    ok = np.isfinite(v)
                    rec[f'class_{k}'] = classify(v[ok], rp.depth.to_numpy()[ok], tol) if ok.any() else None
                    rec[f'obs_matched_{k}'] = float(rp.temperature.to_numpy()[ok].mean()) if ok.any() else np.nan
                prof_rows.append(rec)
        ent_rows.append(row)
    E = pd.DataFrame(ent_rows)
    E = E.apply(verdicts, axis=1)
    return E, pd.DataFrame(prof_rows), pd.DataFrame(series_rows), pd.DataFrame(glamos_rows)


def verdicts(r):
    """Per-cause verdicts by the rules in the header. All three causes warm the ice, so a miss
    where the model is colder than measured refutes all three there."""
    if r.miss in ('near-PMP flip', 'clear miss') and r.err_today < 0:
        msg = f'refuted: model colder than measured ({r.err_today:+.2f} C); this cause would warm it'
        r['cause1_snow'] = r['cause2_water'] = r['cause3_memory'] = msg
        return r
    # (1) snow
    if r.n_mays_band_ice == 0:
        v1 = 'refuted: band firn in all pre-reading Mays (snow inert on firn bands)'
    elif isinstance(r.get('glamos'), str) and pd.notna(r.get('ratio_band')):
        q = r.ratio_band
        v1 = ('supported' if q > SNOW_SUPPORT else 'refuted' if q < SNOW_REFUTE else 'undecidable') + \
             f': model/GLAMOS band winter balance {q:.2f} ({r.glamos_years})'
    else:
        v1 = f'undecidable: no measured winter balance ({r.n_mays_band_ice}/{r.n_mays_total} Mays on ice)'
    # (2) meltwater
    if r.n_years_firn_pre10 == 0:
        v2 = 'refuted: band ice in the 10 years before (no percolation on ice bands)'
    elif r.n_years_wet_to_obs_pre10 > 0:
        v2 = f'possible: water reached {r.d_obs_min:.0f} m in {r.n_years_wet_to_obs_pre10} firn years before (H)'
    else:
        v2 = f'refuted: firn, but water never reached {r.d_obs_min:.0f} m in the 10 years before (H, max {r.wet_max_pre10:.0f} m)'
    # (3) memory: the ice run that ends at the first reading
    lr = r.firn_lost_before_reading
    if not r.reading_months:
        v3 = 'undecidable: no dated reading'
    elif lr == 'firn at reading':
        v3 = 'not applicable: band firn at the first reading'
    elif lr == 'ice since 1972':
        v3 = 'refuted: band ice in every year before the reading'
    elif isinstance(lr, (int, np.integer)):
        yrs = min(int(h) for h in r.reading_months.split(',')) // 100 - lr
        near = r.today > -NEAR_ZERO
        tag = f'firn lost {lr} ({yrs} yr before the first reading), model {r.today:+.2f} vs measured {r.obs:+.2f} C'
        if r.class_obs == r.class_today:
            v3 = 'no miss to explain: ' + tag
        elif yrs <= MEMORY_YEARS and near:
            v3 = 'supported: ' + tag
        elif not near:
            v3 = 'refuted (model already cold): ' + tag
        else:
            v3 = 'undecidable: ' + tag
    else:
        v3 = f'undecidable: {lr}'
    r['cause1_snow'], r['cause2_water'], r['cause3_memory'] = v1, v2, v3
    return r


def raw_readings():
    """Every glenglat reading DataHandler would keep (equilibrium not 'false', from 1950, months
    reaching into the run), with its band, dating and measurement unit (a profile, or a
    calendar month of a logger hole)."""
    b, p, m = load_glenglat()
    p = p[p['eqs'].isin(EQ_KEEP) & (pd.to_datetime(p.date_min, errors='coerce').dt.year >= 1950)]
    dm = {(r.borehole_id, r.id): date_months(r.date_min, r.date_max) for r in p.itertuples()}
    p = p.assign(dm=[dm[(r.borehole_id, r.id)] for r in p.itertuples()])
    p = p[p.dm.notna()]
    p = p.assign(ym0=[d[0] for d in p.dm], ym1=[d[1] for d in p.dm], frac=[d[2] for d in p.dm])
    p = p[(p.ym1 >= WINDOW[0]) & (p.ym0 <= WINDOW[1])]
    n_prof = p.groupby('borehole_id').size()
    p['unit'] = [f'm{r.ym0}' if n_prof[r.borehole_id] > LOGGER_PROFILES else f'p{r.id}' for r in p.itertuples()]
    r = m.merge(p[['borehole_id', 'id', 'eqs', 'date_min', 'date_max', 'ym0', 'ym1', 'frac', 'unit']],
                left_on=['borehole_id', 'profile_id'], right_on=['borehole_id', 'id'])
    r = r.merge(b[['id', 'glacier_name', 'elevation', 'label']].rename(columns={'id': 'borehole_id'}),
                on='borehole_id')
    r = r[r.elevation.notna() & r.temperature.notna() & r.depth.notna()]
    r['band'] = (np.floor(r.elevation / 10) * 10 + 5).astype(int)
    return r


# ---- part B ----------------------------------------------------------------------------------
def load_glenglat():
    b = pd.read_csv(GLENGLAT / 'borehole.csv')
    p = pd.concat([pd.read_csv(GLENGLAT / 'profile.csv'), pd.read_csv(GLENGLAT / 'beer2026/profile.csv')],
                  ignore_index=True)
    m = pd.concat([pd.read_csv(GLENGLAT / 'measurement.csv'),
                   pd.read_csv(GLENGLAT / 'beer2026/measurement.csv')], ignore_index=True)
    p['eqs'] = p.equilibrium.astype(object).where(p.equilibrium.notna(), '').astype(str).str.lower().replace({'nan': '', '<na>': '', 'none': ''})
    return b, p, m


def part_b(gl, readings, bdf, lookup):
    b, p, m = load_glenglat()
    alps = b[b.latitude.between(*ALPS['lat']) & b.longitude.between(*ALPS['lon'])]
    holes = alps[alps.to_bed == True]  # noqa: E712  (column holds True/False/NaN)
    look = pd.read_csv(GLIMS_LOOKUP)
    look['gid'] = look.rgi7_glogem.str.split('.').str[-1]
    glims = dict(zip(look.glims_id, look.gid))
    by_name = {g.glacier_name: g for g in gl}
    run_ids = {f.stem.split('_')[-1] for f in run_dir('today').glob('firnice_snowmay_*.dat')}
    deepest = m.loc[m.groupby(['borehole_id', 'profile_id']).depth.idxmax()]
    deepest = deepest.merge(p[['borehole_id', 'id', 'eqs', 'date_min', 'date_max']],
                            left_on=['borehole_id', 'profile_id'], right_on=['borehole_id', 'id'])
    rows, units = [], []
    for _, h in holes.iterrows():
        bid = int(h.id)
        tol = borehole_tolerance(str(bid), bdf)
        ent = f'{h.glacier_name}@{int(np.floor(h.elevation / 10) * 10 + 5)}'
        ent = ent if ent in by_name else None
        gid = by_name[ent].glacier_id if ent in by_name else glims.get(h.glims_id)
        band_h = int(np.floor(h.elevation / 10) * 10 + 5)
        rec = dict(borehole_id=bid, glacier=h.glacier_name, label=h.label, source=h.source_id,
                   elevation=h.elevation, band_hole=band_h, hole_depth=h.depth, tolerance=tol,
                   entity=ent, rgi=gid, mapping='entity' if ent else ('GLIMS lookup' if gid else 'none'))
        dp = deepest[deepest.borehole_id == bid]
        rec['n_profiles'] = len(dp)
        rec['n_profiles_false'] = int((dp['eqs'] == 'false').sum())
        dp = dp[dp['eqs'].isin(EQ_KEEP)]
        rec['n_profiles_used'] = len(dp)
        if dp.empty:
            rec['status'] = 'no equilibrated profile'
            rows.append(rec)
            continue
        # measurement units: each profile, loggers pooled per calendar month at their modal depth
        dp = dp.assign(dm=[date_months(a, c) for a, c in zip(dp.date_min, dp.date_max)])
        dp = dp[dp.dm.notna()]
        logger = len(dp) > LOGGER_PROFILES
        if logger:
            zmode = dp.depth.mode().iloc[0]
            dp = dp[dp.depth == zmode]
            dp = dp.assign(key=[d[0] for d in dp.dm])
            u = dp.groupby('key').agg(T=('temperature', 'mean'), depth=('depth', 'first'), n=('temperature', 'size'),
                                      date=('date_min', 'min')).reset_index()
            u['dm'] = [(int(k), int(k), 0.5) for k in u.key]
        else:
            u = dp.rename(columns={'temperature': 'T'})[['T', 'depth', 'date_min', 'dm']].rename(columns={'date_min': 'date'})
            u['n'] = 1
        u['pmp'] = pressure_melting_point(u.depth)
        u['frozen_obs'] = u['T'] < u.pmp - tol
        rec['logger'] = logger
        rec['deepest_reading'] = float(u.depth.max())
        rec['gap_to_hole_depth'] = float(h.depth - u.depth.max()) if pd.notna(h.depth) else np.nan
        rec['gap_flag'] = bool(rec['gap_to_hole_depth'] > BED_GAP_FLAG) if pd.notna(h.depth) else None
        rec['T_bed_obs'] = float(u['T'].mean())
        rec['date_first'], rec['date_last'] = str(u.date.min()), str(u.date.max())
        rec['frozen_obs'] = bool((u['T'].mean()) < pressure_melting_point(u.depth.mean()) - tol)
        rec['n_units'] = len(u)
        # model
        if gid is None or gid not in run_ids:
            rec['status'] = 'glacier not run' if gid else 'no RGI7 outline (GLIMS id not in lookup)'
            rows.append(rec)
            continue
        snow_bands = read_snowmay(gid).index.to_numpy()
        g = SimpleNamespace(glacier_id=gid, elevation_band=band_h, glacier_name=f'{h.glacier_name}#{bid}')
        band_m, src = model_band(g, lookup, snow_bands)
        rec['band_model'], rec['model_source'] = band_m, src
        bt = band_table(gid)
        thick = float(bt.thick.get(band_m, np.nan))
        rec['band_thickness'] = thick
        rec['own_bed_node'] = float(NODES[np.argmin(np.abs(NODES - thick))]) if np.isfinite(thick) else np.nan
        u_in = [(i, r) for i, r in u.iterrows() if r.dm[1] >= WINDOW[0] and r.dm[0] <= WINDOW[1]]
        rec['n_units_in_run'] = len(u_in)
        if u_in:
            uin = pd.DataFrame([r for _, r in u_in])
            rec['T_bed_obs_in_run'] = float(uin['T'].mean())
            rec['frozen_obs_in_run'] = bool(uin['T'].mean() < pressure_melting_point(uin.depth.mean()) - tol)
        if not src.startswith('temp_ID'):
            rec['status'] = 'band has no profile output (bedrock summary only)'
            for k in RUNS:
                e, yrs, data = read_summary_file(run_dir(k) / f'temp_bedrock_{gid}.dat')
                ib = int(np.flatnonzero(e == band_m)[0])
                vals = [np.nanmean([data[ib, list(yrs).index(hyd_year(ym))] for ym in
                                    reading_months(*r.dm)[0] if hyd_year(ym) in list(yrs)]) for _, r in u_in]
                rec[f'Tbed_{k}'] = float(np.nanmean(vals)) if vals else np.nan
                rec[f'frozen_{k}'] = (bool(rec[f'Tbed_{k}'] < BEDROCK_CAP - 1e-3)
                                      if np.isfinite(rec[f'Tbed_{k}']) else None)
            rec['col_depth'] = rec['own_bed_node']
            rows.append(rec)
            continue
        if not u_in:
            rec['status'] = 'no profile inside the run window'
        else:
            rec['status'] = 'compared'
        for k in RUNS:
            _, _, df = profile(k, gid, src)
            zb, _ = bottom_node(df)
            rec['col_depth'] = zb
            if not u_in:
                continue
            vb, vz = [], []
            for _, r in u_in:
                mo = [reading_months(*r.dm)]
                vb.append(model_at(df, [zb], mo)[0])
                vz.append(model_at(df, [r.depth], mo)[0])
            rec[f'Tbed_{k}'] = float(np.nanmean(vb))
            rec[f'T_at_reading_depth_{k}'] = float(np.nanmean(vz)) if np.isfinite(vz).any() else np.nan
            rec[f'frozen_{k}'] = bool(rec[f'Tbed_{k}'] < pressure_melting_point(zb) - tol)
            zr = float(np.mean([r.depth for _, r in u_in]))
            ta = rec[f'T_at_reading_depth_{k}']
            rec[f'frozen_atdepth_{k}'] = bool(ta < pressure_melting_point(zr) - tol) if np.isfinite(ta) else None
            if k == 'today':
                for (i, r), a, z in zip(u_in, vb, vz):
                    units.append(dict(borehole_id=bid, glacier=h.glacier_name, date=r.date, depth=r.depth,
                                      T_obs=r['T'], frozen_obs=bool(r.frozen_obs), n=int(r.n),
                                      col_depth=zb, Tbed_today=a, T_at_depth_today=z))
        rec['deepened'] = bool(rec['col_depth'] > rec['own_bed_node'] + 0.5)
        rows.append(rec)
    return pd.DataFrame(rows), pd.DataFrame(units)


def confusion(B, col='frozen'):
    """col='frozen': model bottom node; 'frozen_atdepth': model at the measured deepest reading."""
    c = B[B.status == 'compared']
    out = {}
    for k in RUNS:
        cc = c[c[f'{col}_{k}'].notna()]
        t = pd.crosstab(cc.frozen_obs_in_run.map({True: 'obs frozen', False: 'obs temperate'}),
                        cc[f'{col}_{k}'].astype(bool).map({True: 'model frozen', False: 'model temperate'}))
        out[k] = {f'{i} / {j}': int(t.loc[i, j]) if i in t.index and j in t.columns else 0
                  for i in ('obs frozen', 'obs temperate') for j in ('model frozen', 'model temperate')}
    return out


# ---- part C ----------------------------------------------------------------------------------
def part_c(lookup):
    """temp_bedrock = annual max of tl_fit[nfit], the below-bed fill slot, clipped by the universal
    PMP clamp at the PMP of the 519 m slot; checked against the temp_ID bottom node."""
    run_ids = sorted(f.stem.split('_')[-1] for f in run_dir('today').glob('firnice_snowmay_*.dat'))
    lookup_depth = {}
    for ln in (REPO / 'data/bayescal' / RUNS['today'] / 'training_runs/glenglat_borehole_elevations.dat').read_text().splitlines():
        a = ln.split()
        if a and not a[0].startswith('#'):
            lookup_depth.setdefault(a[0], [t.split(':')[1] if ':' in t else -1 for t in a[1:]])
    checks, types = [], []
    for gid in run_ids:
        e, yrs, data = read_summary_file(run_dir('today') / f'temp_bedrock_{gid}.dat')
        yrs = list(yrs)
        for pf in sorted(run_dir('today').glob(f'temp_ID*_{gid}.dat')):
            elev, _, df = profile('today', gid, re.match(r'(temp_ID\d+)_', pf.name).group(1))
            zb, Tb = bottom_node(df)
            if Tb is None or elev not in set(e):
                continue
            ib = int(np.flatnonzero(e == elev)[0])
            amax = pd.Series(Tb, index=df.year.to_numpy()).groupby(level=0).max()
            pred = np.minimum(amax.to_numpy(), BEDROCK_CAP)
            got = np.array([data[ib, yrs.index(y)] for y in amax.index])
            k = int(re.match(r'temp_ID(\d+)_', pf.name).group(1)) - 1
            md = float(lookup_depth[gid][k]) if gid in lookup_depth and k < len(lookup_depth[gid]) else -1.0
            th = float(band_table(gid).thick.get(int(elev), np.nan))
            own = float(NODES[np.argmin(np.abs(NODES - th))]) if np.isfinite(th) else np.nan
            expect = own if md <= th else float(NODES[NODES >= md][0]) if (NODES >= md).any() else float(NODES[-1])
            checks.append(dict(rgi=gid, band=int(elev), col_depth=zb, lookup_depth=md, band_thickness=th,
                               own_node=own, expected_col=expect, col_ok=bool(np.isclose(expect, zb)),
                               max_abs_diff=float(np.nanmax(np.abs(pred - got))),
                               n_years=len(got), share_at_cap=float(np.mean(np.isclose(got, BEDROCK_CAP, atol=6e-4)))))
        bt = band_table(gid)
        for y in TYPE_YEARS:
            if y not in yrs:
                continue
            col = data[:, yrs.index(y)]
            ok = np.isfinite(col)
            el, T = e[ok], col[ok]
            area = bt.area.reindex(el).fillna(0).to_numpy()
            frozen = T < BEDROCK_CAP - 1e-3
            o = np.argsort(el)
            fr = frozen[o]
            if fr.all():
                t = 'cold bed'
            elif not fr.any():
                t = 'temperate bed'
            elif fr[0] and not fr[-1]:
                t = 'polythermal: frozen front'
            elif fr[-1] and not fr[0]:
                t = 'polythermal: frozen upper part'
            else:
                t = 'mixed'
            types.append(dict(rgi=gid, year=y, type=t, frozen_area_share=float(area[frozen].sum() / area.sum()),
                              n_bands=int(ok.sum()), n_frozen=int(frozen.sum()),
                              frozen_bands=compress_ranges(el[frozen])))
    return pd.DataFrame(checks), pd.DataFrame(types)


def compress_ranges(el):
    el = sorted(int(x) for x in el)
    if not el:
        return ''
    out, a, b = [], el[0], el[0]
    for x in el[1:]:
        if x == b + 10:
            b = x
        else:
            out.append(f'{a}-{b}' if a != b else f'{a}')
            a = b = x
    out.append(f'{a}-{b}' if a != b else f'{a}')
    return ','.join(out)


# ---- main ------------------------------------------------------------------------------------
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P = pd.read_csv(POINTS, dtype={'rgi': str})
    d = pickle.load(open(CACHE, 'rb'))
    gl, readings = d['glaciers'], d['readings']
    bdf = load_borehole_table(GLENGLAT)
    lookup = read_lookup()
    log = []
    say = lambda s='': (print(s), log.append(s))

    # the mass-balance side is identical across runs (firn line, snow)
    same = all((run_dir(k) / f.name).read_bytes() == f.read_bytes()
               for k in RUNS for f in run_dir('today').glob('firnice_firnline_*.dat'))
    say(f'firn-line files identical across runs: {same}')

    raw = raw_readings()
    cnt = []
    for g in gl:
        rr = raw[(raw.glacier_name == g.base_glacier_name) & (raw.band == int(g.elevation_band))]
        cnt.append(dict(entity=g.glacier_name, rgi=g.glacier_id, raw=len(rr), cache=int(np.sum(g.n_readings)),
                        raw_profiles=rr.groupby(['borehole_id', 'id']).ngroups,
                        raw_deep=int((rr.depth > DEEP).sum())))
    cnt = pd.DataFrame(cnt).assign(missing=lambda x: x.raw - x.cache)
    cnt.to_csv(OUT / 'entity_reading_counts.csv', index=False)
    say(f'readings DataHandler would keep vs the cache: {int((cnt.missing > 0).sum())} of {len(cnt)} entities '
        f'miss {int(cnt.missing.clip(lower=0).sum())} of {int(cnt.raw.sum())} readings; '
        f'{int((cnt.missing < 0).sum())} entities have more in the cache')
    E, PR, S, G = part_a(P, gl, raw, bdf, lookup)
    E.to_csv(OUT / 'zone_entities.csv', index=False, float_format='%.3f')
    PR.to_csv(OUT / 'zone_profiles.csv', index=False, float_format='%.3f')
    S.to_csv(OUT / 'zone_T15_T20_by_year.csv', index=False, float_format='%.3f')
    G.to_csv(OUT / 'zone_glamos_winter.csv', index=False, float_format='%.3f')
    deep = E[E.n_deep > 0]
    say(f'\nA. zone {ZONE}: {len(E)} entities, {len(deep)} with points below {DEEP:.0f} m; '
        f'regime right (today) {int((deep.miss == "hit").sum())}/{len(deep)}')
    cols = ['entity', 'band_model', 'n_deep', 'obs', 'today', 'F', 'H', 'FH', 'class_obs', 'class_today',
            'class_F', 'class_H', 'class_FH', 'miss']
    say(deep[cols].to_string(index=False, float_format='%.2f'))
    cols = ['entity', 'snow_may_pre5', 'n_mays_band_ice', 'n_mays_total', 'n_years_firn_pre10',
            'n_years_wet_to_obs_pre10', 'firn_loss_year', 'firn_at_reading', 'T20_loss-1', 'T20_loss+5',
            'T20_loss+10', 'T20_at_reading']
    say(E[[c for c in cols if c in E.columns]].to_string(index=False, float_format='%.2f'))
    for _, r in E.iterrows():
        say(f'{r.entity} [{r.miss}]\n   1 {r.cause1_snow}\n   2 {r.cause2_water}\n   3 {r.cause3_memory}')
    if 'ratio_band' in E:
        say(E[E.ratio_band.notna()][['entity', 'glamos', 'glamos_years', 'glamos_overlap_pre5',
                                     'model_band_mean', 'glamos_band_mean', 'ratio_band',
                                     'ratio_glacier']].to_string(index=False, float_format='%.2f'))
    say('\nprofiles of the zone (deep points), measured vs model at the readings\' dates:')
    say(PR[['entity', 'borehole_id', 'unit', 'date', 'n', 'depth_min', 'depth_max', 'obs', 'today', 'F', 'H', 'FH',
            'class_obs', 'class_today']].to_string(index=False, float_format='%.2f'))

    B, U = part_b(gl, readings, bdf, lookup)
    B.to_csv(OUT / 'bed_holes.csv', index=False, float_format='%.3f')
    U.to_csv(OUT / 'bed_units_today.csv', index=False, float_format='%.3f')
    conf = dict(bottom_node=confusion(B), at_reading_depth=confusion(B, 'frozen_atdepth'))
    say(f'\nB. {len(B)} Alpine holes with to_bed, {B.glacier.nunique()} glaciers; status:')
    say(B.status.value_counts().to_string())
    cols = ['borehole_id', 'glacier', 'label', 'band_hole', 'band_model', 'hole_depth', 'deepest_reading',
            'band_thickness', 'col_depth', 'deepened', 'T_bed_obs_in_run', 'frozen_obs_in_run', 'Tbed_today',
            'frozen_today', 'frozen_F', 'frozen_H', 'frozen_FH', 'status']
    say(B[[c for c in cols if c in B.columns]].to_string(index=False, float_format='%.2f'))
    say(json.dumps(conf, indent=1))

    ela = []
    for gid, name in GLAMOS_MAP.items():
        gb = glamos_bands(name)
        fl = read_firnline('today', gid)
        for r in gb['rows']:
            if r['year'] in fl.index:
                ela.append(dict(rgi=gid, glamos=name, year=r['year'], ela_glamos=r['ela'], glamos_hmax=r['hmax'],
                                ela_model=float(fl.loc[r['year'], 'ela']),
                                firnline_model=float(fl.loc[r['year'], 'firnline'])))
    ela = pd.DataFrame(ela)
    ela.to_csv(OUT / 'glamos_ela.csv', index=False)
    say('\nELA, GLAMOS against the model (annual ELA and the 5-year firn line), mean over common years:')
    say(ela.groupby('glamos').agg(years=('year', lambda y: f'{y.min()}-{y.max()} ({len(y)})'),
                                  ela_glamos=('ela_glamos', 'mean'), ela_model=('ela_model', 'mean'),
                                  firnline_model=('firnline_model', 'mean')).round(0).to_string())
    snow = G.groupby('glamos').agg(model_glacier=('model_glacier_wide', 'mean'),
                                   glamos_glacier=('glamos_bw_glacier', 'mean'))
    say('glacier-wide winter balance, model 1 May snow / GLAMOS Bw [m w.e.]:')
    say(snow.round(2).to_string())

    C, T = part_c(lookup)
    C.to_csv(OUT / 'bedrock_file_check.csv', index=False, float_format='%.4f')
    T.to_csv(OUT / 'glacier_bed_type.csv', index=False, float_format='%.3f')
    say(f'\nC. temp_bedrock vs min(annual max of temp_ID bottom node, {BEDROCK_CAP}): '
        f'max |diff| {C.max_abs_diff.max():.4f} C over {len(C)} profile bands; column depth as '
        f'expected from band thickness and lookup depth at {int(C.col_ok.sum())}/{len(C)}; '
        f'deepened beyond the band thickness at {int((C.expected_col > C.own_node + 0.5).sum())}')
    say(C[~C.col_ok].to_string(index=False))
    say(T.pivot(index='rgi', columns='year', values='type').to_string())

    meta = dict(selections=dict(zone=ZONE, deep=DEEP, near_pmp=NEAR_PMP, snow_years=SNOW_YEARS,
                                wet_years=WET_YEARS, wet_run=WET_RUN, glamos=GLAMOS_MAP,
                                snow_support=SNOW_SUPPORT, snow_refute=SNOW_REFUTE, near_zero=NEAR_ZERO,
                                memory_years=MEMORY_YEARS, alps=ALPS, eq_keep=sorted(EQ_KEEP),
                                window=WINDOW, bed_gap_flag=BED_GAP_FLAG, logger_profiles=LOGGER_PROFILES,
                                bedrock_cap=BEDROCK_CAP, default_tolerance=DEFAULT_TOLERANCE),
                runs=RUNS, points=str(POINTS), cache=str(CACHE), confusion=conf,
                firnline_identical=same)
    (OUT / 'summary.json').write_text(json.dumps(meta, indent=1, default=str))
    (OUT / 'zone3000_bed_check.log').write_text('\n'.join(log) + '\n')


if __name__ == '__main__':
    main()
