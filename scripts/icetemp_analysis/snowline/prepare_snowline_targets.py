"""Prepare transient-snowline (TSLA) observation files and a per-glacier target table
for the CentralEurope_fair calibration glaciers (snowline calibration after von der Esch et al.).

Inputs
  Sommer et al. (ESSD preprint essd-2026-35) Alpine SLA product, RGI7 IDs
  GloGEM RGI7 band files (area, hypsometry), Hugonnet 2000-2020 RGI7 geodetic balances,
  GLAMOS *_fix.dat glacier-wide balances, GLAMOS<->RGI7 matching in targets/glamos_vs_geodetic.csv

Outputs (in OUT)
  tsla_csv/tsla_RGI70-11.xxxxx.csv   one clean observation file per glacier
  snowline_2000_2020/11_snowline_median.dat   all glaciers, read by GloGEM's snowline calibration
  eos_sla.csv                        Sommer end-of-summer SLA per glacier and year, with SO/TO codes
  tsla_counts.csv                    observation counts per glacier and QA level
  targets.csv                        target table for all catchment glaciers

Usage: python prepare_snowline_targets.py [--qa-min 0.6666]
"""
import argparse
import os

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA = f'{REPO}/data/snowlines'
OUT = DATA

SOM = f'{DATA}/sommer/transient_snow_line_altitudes_of_glaciers_europ_alps'
MATCH = f'{DATA}/reference/glamos_vs_geodetic.csv'
CATCH = '/itet-stor/jabeer/glogem/data/catchments/RGI11_CentralEurope_fair.dat'
BANDS = '/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope'
GEOD = '/itet-stor/jabeer/glogem/data/geodetic/RGIv7.0/aggregated_2000_2020/11_mb_glspec.dat'
GLAMOS = '/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier'

Y0, Y1 = 2000, 2019          # calibration window, hydrological years (Oct Y-1 .. Sep Y)
OBS_END = '2025-12-31'
EDGE_TOL = 20.0              # m, tolerance for the inferred full-snow flag
MIN_TSLA_AREA = 1.0          # km2, Sommer product only covers glaciers >= 1 km2
QA_LEVELS = {'qa_ge0': 0.0, 'qa_ge4of6': 4 / 6 - 1e-6, 'qa_ge5of6': 5 / 6 - 1e-6, 'qa_eq1': 1 - 1e-6}

NAMES = {
    '00245': 'Pelvoux', '00773': 'Bossons', '00774': 'Taconnaz', '00777': 'Tete Rousse',
    '00786': 'Mont Blanc summit', '01039': 'Tortin', '01112': 'Hohsaas', '01127': 'Gruben',
    '01150': 'Alphubel N', '01195': 'Colle Gnifetti area', '01205': 'Lys', '01225': 'Gorner/Grenz',
    '01230': 'Monte Rosa E face', '01265': 'Altels', '01309': 'Sex Rouge', '01565': 'Sphinxgrat',
    '02103': 'Corvatsch', '02216': 'Morteratsch', '02596': 'Gr. Aletsch', '02611': 'Fiescher',
    '02669': 'Titlis', '02821': 'Ortles (Alta)', '03116': 'Hintereis',
}


def catchment_ids():
    ids = [ln.strip() for ln in open(CATCH) if ln.strip().startswith('RGI70-11.')]
    return [i[-5:] for i in ids]


def read_bands(gid):
    """RGI7 area from the header and band hypsometry (2003 inventory geometry)."""
    with open(f'{BANDS}/{gid}.dat') as f:
        hdr = [next(f) for _ in range(3)]
    area = float(hdr[2].split(':')[1].split()[0])
    b = np.loadtxt(f'{BANDS}/{gid}.dat', skiprows=5, usecols=(1, 2, 3))
    return area, b[:, 0], b[:, 1], b[:, 2]


def hyps_quantile(lo, hi, a, q):
    """Elevation below which a fraction q of the glacier area lies."""
    ca = np.concatenate([[0], np.cumsum(a) / a.sum()])
    return np.interp(q, ca, np.concatenate([[lo[0]], hi]))


def read_eos():
    v = pd.read_csv(f'{SOM}/02_end_of_summer_glacier_snow_lines/Alps_Glacier_EoS_SLA_2000-2025.csv', dtype=str)
    d = pd.read_csv(f'{SOM}/02_end_of_summer_glacier_snow_lines/Alps_Glacier_EoS_SLA_DOY_2000-2025.csv', dtype=str)
    return v.set_index(v.Year.astype(int)), d.set_index(d.Year.astype(int))


def hydro_year(dates):
    return dates.dt.year + (dates.dt.month >= 10).astype(int)


def prepare_tsla(gid, qa_min, eos_v, eos_d):
    fn = f'{SOM}/01_glacier_snow_lines/11-01_RGI2000-v7.0-G-11-{gid}.csv'
    if not os.path.exists(fn):
        return None, None, None
    d = pd.read_csv(fn, parse_dates=['Date'])
    d = d[(d.Date >= '2000-01-01') & (d.Date <= OBS_END)].sort_values('Date').reset_index(drop=True)
    d['hydro_year'] = hydro_year(d.Date)
    cal = d.hydro_year.between(Y0, Y1)
    cnt = {'gid': gid}
    for k, q in QA_LEVELS.items():
        cnt[f'n_{k}'] = int((d.QA_FLAG >= q).sum())
        cnt[f'n_{k}_cal'] = int(((d.QA_FLAG >= q) & cal).sum())
    qa_ok = (d.QA_FLAG >= qa_min) & cal
    per_yr = d[qa_ok].groupby('hydro_year').size().reindex(range(Y0, Y1 + 1), fill_value=0)
    cnt['obs_per_yr_cal_median'] = float(per_yr.median())
    cnt['obs_per_yr_2000_2012'] = float(per_yr.loc[2000:2012].mean())
    cnt['obs_per_yr_2013_2019'] = float(per_yr.loc[2013:2019].mean())
    cnt['yrs_without_obs_cal'] = int((per_yr == 0).sum())
    cnt['sigma_median'] = float(d.SLA_ERROR[qa_ok].median())
    cnt['sigma_p90'] = float(d.SLA_ERROR[qa_ok].quantile(0.9))
    cnt['dlt_elev_median'] = float((d.SLA_CORR - d.SLA)[qa_ok].median())

    d = d[d.QA_FLAG >= qa_min].copy()
    # inferred edge cases: SLA is the 10th percentile of snow-pixel elevations, so a fully
    # snow-covered glacier gives SLA ~ hypsometric 10th percentile; SLA above the 95th
    # percentile means snow on <~5 % of the area
    _, lo, hi, a = read_bands(gid)
    z10, z95 = hyps_quantile(lo, hi, a, 0.10), hyps_quantile(lo, hi, a, 0.95)
    flag = np.where(d.SLA <= z10 + EDGE_TOL, 'full_snow_inferred',
                    np.where(d.SLA >= z95, 'snow_free_inferred', ''))
    # Sommer EoS codes: SO = statistical outlier (+-2 sd), TO = threshold outlier (+-400 m)
    col = f'11-{gid}'
    win = (d.Date.dt.strftime('%m-%d') >= '07-15') & (d.Date.dt.strftime('%m-%d') <= '09-30')
    eos_code = np.array([eos_v.loc[y, col] if (w and col in eos_v and eos_v.loc[y, col] in ('SO', 'TO')) else ''
                         for y, w in zip(d.Date.dt.year, win)], dtype=object)
    eos_flag = np.where(eos_code != '', 'eos_' + eos_code.astype(str), '')
    flag = [';'.join(s for s in (f, e) if s) for f, e in zip(flag, eos_flag)]
    cnt['n_full_snow_inferred'] = int(sum('full_snow' in f for f in flag))
    cnt['n_snow_free_inferred'] = int(sum('snow_free' in f for f in flag))
    cnt['z_hyps_p10'], cnt['z_hyps_p95'] = round(z10, 1), round(z95, 1)
    out = pd.DataFrame({
        'date': d.Date.dt.strftime('%Y-%m-%d'), 'hydro_year': d.hydro_year, 'month': d.Date.dt.month,
        'sla_m': d.SLA_CORR.round(2), 'sigma_m': d.SLA_ERROR.round(2), 'qa': d.QA_FLAG.round(4),
        'sensor': d.Sensor, 'flag': flag, 'sla_raw_m': d.SLA.round(2)})
    eos = []
    if col in eos_v:
        for y in eos_v.index:
            v, doy = eos_v.loc[y, col], eos_d.loc[y, col]
            code = v if v in ('SO', 'TO') else ('' if pd.notna(v) else 'NA')
            eos.append(dict(rgi7=f'RGI70-11.{gid}', year=y, eos_sla_m=pd.to_numeric(v, errors='coerce'),
                            eos_doy=pd.to_numeric(doy, errors='coerce'), code=code))
    return out, cnt, eos


def doy365(dates):
    """Day of year on GloGEM's daily calendar (365 days; 29 Feb taken as 28 Feb)."""
    return dates.dt.dayofyear - (dates.dt.is_leap_year & (dates.dt.dayofyear >= 60)).astype(int)


def write_glogem(tsla, qa_min):
    """All glaciers in one file, in the format of GloGEM's read_snowlinedata.pro (A. von der Esch):
    one header line, then: RGI id, snow line [m a.s.l.], stdev [m], year, day of year (365-day calendar).
    11_snowline_median.dat holds SLA_CORR (on the surface of the observation date), _raw.dat the SLA on
    the SRTM DEM of February 2000, closer to GloGEM's fixed inventory geometry."""
    os.makedirs(f'{OUT}/snowline_2000_2020', exist_ok=True)
    n = sum(len(df) for df in tsla.values())
    for suffix, col, label in (('', 'sla_m', 'SLA_CORR'), ('_raw', 'sla_raw_m', 'SLA (raw, SRTM 2000)')):
        with open(f'{OUT}/snowline_2000_2020/11_snowline_median{suffix}.dat', 'w') as f:
            f.write(f'RGI-ID SLA_m sigma_m year doy365  (Sommer et al. TSLA, {label} / SLA_ERROR, QA>={qa_min:.3f}, '
                    f'2000-2025, {len(tsla)} glaciers, {n} observations)\n')
            for gid, df in tsla.items():
                d = pd.to_datetime(df.date)
                for s, e, y, j in zip(df[col], df.sigma_m, d.dt.year, doy365(d)):
                    f.write(f'RGI70-11.{gid} {s:8.1f} {e:6.1f} {y:4d} {j:3d}\n')


def read_fix(short):
    rows = []
    with open(f'{GLAMOS}/{short}_fix.dat', encoding='latin-1') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            a = line.split()
            f = lambda s: np.nan if s.lower() == 'nan' else float(s.rstrip('.')) / 1000
            rows.append(dict(meth=int(a[0]), yr=int(a[4][:4]), Bw=f(a[5]), Ba=f(a[6])))
    return pd.DataFrame(rows)


def glamos_period_adjusted():
    """Mean Ba/Bw over each glacier's 2000-2019 years, minus the regional anomaly of those
    years (from GLAMOS series complete over 2000-2019), as in targets/compare.py."""
    import glob
    ser = {os.path.basename(f)[:-8]: read_fix(os.path.basename(f)[:-8]) for f in sorted(glob.glob(f'{GLAMOS}/*_fix.dat'))}
    reg = {}
    for v in ('Ba', 'Bw'):
        full = {s: d.dropna(subset=[v]).set_index('yr')[v] for s, d in ser.items()}
        full = {s: x.loc[Y0:Y1] for s, x in full.items() if set(range(Y0, Y1 + 1)) <= set(x.index)}
        ref = pd.DataFrame(full)
        reg[v] = (ref - ref.mean()).mean(axis=1)
        reg[v + '_n'] = ref.shape[1]
    return ser, reg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qa-min', type=float, default=4 / 6 - 1e-6,
                    help='QA threshold; QA values are k/6, so 0.67 would drop the 4/6 level')
    qa_min = ap.parse_args().qa_min
    os.makedirs(f'{OUT}/tsla_csv', exist_ok=True)
    ids = catchment_ids()
    eos_v, eos_d = read_eos()

    counts, eos_all, ntsla, tsla = [], [], {}, {}
    for gid in ids:
        df, cnt, eos = prepare_tsla(gid, qa_min, eos_v, eos_d)
        if df is None:
            ntsla[gid] = (0, 0)
            continue
        df.to_csv(f'{OUT}/tsla_csv/tsla_RGI70-11.{gid}.csv', index=False)
        tsla[gid] = df
        counts.append(cnt)
        eos_all += eos
        ntsla[gid] = (len(df), int(df.hydro_year.between(Y0, Y1).sum()))
    write_glogem(tsla, qa_min)
    counts = pd.DataFrame(counts)
    counts.insert(1, 'name', counts.gid.map(NAMES))
    counts.to_csv(f'{OUT}/tsla_counts.csv', index=False, float_format='%.1f')
    pd.DataFrame(eos_all).to_csv(f'{OUT}/eos_sla.csv', index=False)

    # geodetic target 2000-2020 and GloGEM's outlier replacement (read_geodeticdata.pro)
    g = pd.read_csv(GEOD, sep=r'\s+', skiprows=3, header=None, names=['rid', 'area', 'lat', 'lon', 'B', 'errB', 'flag'])
    g['gid'] = g.rid.str[-5:]
    g = g.set_index('gid')
    m = g.flag == 1
    regmean = (g.area[m] * g.B[m]).sum() / g.area[m].sum()
    mu, sd = g.B[m].mean(), g.B[m].std()
    g['B_used'] = np.where(((g.B < mu - 2 * sd) | (g.B > mu + 2 * sd)) & (g.area < 20), regmean, g.B)

    match = pd.read_csv(MATCH, dtype={'rgi7': str})
    match['rgi7'] = match.rgi7.str.zfill(5)
    ser, reg = glamos_period_adjusted()

    rows = []
    for gid in ids:
        area = read_bands(gid)[0]
        rec = dict(rgi7=f'RGI70-11.{gid}', name=NAMES.get(gid, ''), area_rgi7_km2=area,
                   geod_B=g.B.get(gid), geod_err=g.errB.get(gid), geod_flag=int(g.flag.get(gid)),
                   geod_B_used=g.B_used.get(gid),
                   n_tsla=ntsla[gid][0], n_tsla_cal=ntsla[gid][1])
        mm = match[match.rgi7 == gid]
        if len(mm):
            r0 = mm.iloc[0]
            d = ser[r0.short]
            ov = d[d.yr.between(Y0, Y1)]
            oa, ow = ov.dropna(subset=['Ba']), ov.dropna(subset=['Bw'])
            rec.update(glamos=r0.short, glamos_years=f'{oa.yr.min()}-{oa.yr.max()}', glamos_n_Ba=len(oa),
                       glamos_Ba=oa.Ba.mean() - reg['Ba'].loc[oa.yr].mean(),
                       glamos_Ba_raw=oa.Ba.mean(),
                       glamos_n_Bw=len(ow),
                       glamos_n_Bw_measured=int(ow.meth.isin([1, 2, 4, 5]).sum()),
                       glamos_Bw=(ow.Bw.mean() - reg['Bw'].loc[ow.yr].mean()) if len(ow) else np.nan,
                       glamos_meth=';'.join(f'{k}x{v}' for k, v in sorted(oa.meth.value_counts().items())),
                       glamos_n_seasonal_stakes=int(oa.meth.isin([1, 4]).sum()),
                       area_glamos_km2=r0.A_glamos_grid, area_ratio=r0.A_glamos_grid / area)
            rec['outline_mismatch'] = bool(rec['area_ratio'] < 0.7 or rec['area_ratio'] > 1.3)
        else:
            rec.update(glamos='', outline_mismatch=False)
        has_tsla = rec['n_tsla_cal'] > 0 and area >= MIN_TSLA_AREA
        has_gl = bool(rec['glamos'])
        if has_tsla:
            c = 'geodetic + TSLA'
            if has_gl:
                c += ('; GLAMOS partial outline (ratio %.2f), not used' % rec['area_ratio'] if rec['outline_mismatch']
                      else '; GLAMOS for validation')
        elif has_gl and not rec['outline_mismatch']:
            c = 'GLAMOS Ba' + (' + Bw' if rec['glamos_n_Bw_measured'] else ' (Bw not measured)') + ' instead of geodetic'
            if rec['glamos_n_seasonal_stakes'] + int(oa.meth.isin([3, 6]).sum()) == 0:
                c += '; GLAMOS method unknown (0)'
        elif has_gl:
            c = 'geodetic only; GLAMOS partial outline (ratio %.2f) flagged' % rec['area_ratio']
        else:
            c = 'geodetic only'
        rec['constraint'] = c
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.to_csv(f'{OUT}/targets.csv', index=False, float_format='%.3f')

    pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_colwidth', 70)
    print(f'QA threshold {qa_min:.4f}; GLAMOS complete-series reference: Ba n={reg["Ba_n"]}, Bw n={reg["Bw_n"]}')
    print(counts.to_string(index=False))
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
