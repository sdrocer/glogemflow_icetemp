"""Aggregate firn-model output to GloGEM RGI7 10 m bands and compare with the thermal model."""
import glob, os, re, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = os.environ.get('FIRNMODEL_DIR', '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel')   # port_out/ in, band_tables_v2/ out
from firn_change_port import read_agr

MB = '/scratch_net/vierzack06_fourth/mhuss/MassBalance_all/'
BANDS = '/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope/'
TH = ('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/bayescal/'
      'fair_FRH_r0.85_i1.017_a1.00_w0/training_runs/monthly/CentralEurope/PAST/firnice_temperature/')
RHO_FIRN = 600.
MINCELLS = 3
GL = {'corvatsch': '02103', 'hohsaas': '01112', 'aletsch': '02596', 'fiescher': '02611',
      'gorner': '01225', 'morteratsch': '02216', 'sexrouge': '01309', 'tortin': '01039',
      'alphubel': '01150'}
ZONES = {'corvatsch': (3295, 3365), 'hohsaas': (3275, 3345), 'aletsch': (3415, 3535)}


def band_tables(g, gid):
    r = dict(np.load(f'{S}/port_out/port_{g}.npz'))
    years = r['years']; mask = r['glacier_mask']
    ref, h = read_agr(f'{MB}{g}/massbalance/results/massbalmeas_2025.grid')
    dems = []
    for f in sorted(glob.glob(f'{MB}{g}/topo/dhm_*.grid')):
        z, hz = read_agr(f)
        if z.shape == ref.shape and abs(hz['xll'] - h['xll']) < 1 and abs(hz['yll'] - h['yll']) < 1:
            dems.append((int(re.findall(r'(\d{4})', os.path.basename(f))[-1]), z, os.path.basename(f)))
    cell_km2 = h['cs'] ** 2 / 1e6
    bd = np.loadtxt(f'{BANDS}{gid}.dat', skiprows=5)
    rgi = pd.DataFrame({'band_mid': (bd[:, 1] + bd[:, 2]) / 2, 'rgi_area_km2': bd[:, 3]})
    fl = np.loadtxt(f'{TH}firnice_firnline_{gid}.dat', comments='#')
    fl = pd.DataFrame(fl, columns=['year', 'ela', 'aar', 'firnline_th', 'n_f2i', 'n_i2f'])
    wd = pd.DataFrame(np.loadtxt(f'{TH}firnice_wetdepth_{gid}.dat', comments='#'),
                      columns=['band_mid', 'year', 'wet_depth_m', 'water_lost_mwe'])
    va = np.loadtxt(f'{TH}firnice_temp_validation_{gid}.dat', comments='#')
    zp = pd.DataFrame({'band_mid': va[:, 0], 'z_perm_m': va[:, 9]})
    top_th = va[:, 0].max()
    rows = []
    for i, y in enumerate(years):
        dy, z, dn = min(dems, key=lambda d: (abs(d[0] - y), -d[0]))
        mid = np.floor(z / 10.) * 10 + 5
        ft = r['firntot'][i]; dur = r['dur'][i]
        ftp = np.where(ft > 0, ft, 0.)
        fc = mask & (dur >= 2)
        for bm in np.unique(mid[mask]):
            sel = mask & (mid == bm)
            n = sel.sum(); nf = (sel & fc).sum()
            rows.append(dict(glacier=g, rgi7_id=gid, year=y, band_mid=bm, dem=dn, glamos_cells=n,
                             glamos_area_km2=n * cell_km2, fm_firn_cells=nf, fm_firn_frac=nf / n,
                             fm_firn_mwe_firncells=ft[sel & fc].mean() if nf else np.nan,
                             fm_firn_mwe_allcells=ftp[sel].mean(),
                             mb_mwe=np.nanmean(np.where(r['bal'][i][sel] != -9999, r['bal'][i][sel], np.nan))
                             if r['present'][i] else np.nan,
                             mb_year_present=bool(r['present'][i])))
    fmb = pd.DataFrame(rows)
    fmb['fm_firn_m_firncells'] = fmb.fm_firn_mwe_firncells * 1000. / RHO_FIRN
    # full band x year grid on the RGI7 bands
    grid = pd.MultiIndex.from_product([years, rgi.band_mid], names=['year', 'band_mid']).to_frame(index=False)
    out = grid.merge(rgi, on='band_mid').merge(fmb.drop(columns=['glacier', 'rgi7_id']), on=['year', 'band_mid'], how='left')
    out = out.merge(fl[['year', 'firnline_th', 'ela', 'aar']], on='year', how='left')
    out['th_firn'] = np.where(out.firnline_th.isna(), np.nan, (out.band_mid >= out.firnline_th).astype(float))
    out = out.merge(wd, on=['band_mid', 'year'], how='left').merge(zp, on='band_mid', how='left')
    out.insert(0, 'glacier', g); out.insert(1, 'rgi7_id', gid)
    out['glamos_cells'] = out.glamos_cells.fillna(0).astype(int)
    return out, top_th


def fm_line(d):
    v = d[(d.glamos_cells >= MINCELLS) & (d.fm_firn_frac >= 0.5)]
    return v.band_mid.min() if len(v) else np.nan


def summarise(out, g):
    res = []
    yrs = sorted(out.year.unique())
    ann = []
    for y in yrs:
        d = out[out.year == y]
        ov = d[(d.glamos_cells >= MINCELLS) & (d.rgi_area_km2 > 0)]
        A = ov.rgi_area_km2.sum()
        th = ov.th_firn
        if th.isna().all():
            gg_fm_ice = fm_gg_ice = np.nan
        else:
            gg_fm_ice = (ov.rgi_area_km2 * th * (1 - ov.fm_firn_frac)).sum()
            fm_gg_ice = (ov.rgi_area_km2 * (1 - th) * ov.fm_firn_frac).sum()
        th_line = d.firnline_th.iloc[0]
        ann.append(dict(year=y, fm_line=fm_line(d), th_line=th_line if th_line <= d.band_mid.max() else np.inf,
                        ov_area=A, gg_firn_fm_ice=gg_fm_ice, fm_firn_gg_ice=fm_gg_ice,
                        mb_present=bool(d.mb_year_present.fillna(False).any())))
    ann = pd.DataFrame(ann)
    ann['decade'] = (ann.year // 10) * 10
    for dec, a in ann.groupby('decade'):
        a2 = a[a.mb_present]
        res.append(dict(glacier=g, decade=f'{dec}s', n_mb_years=len(a2),
                        fm_line_mean=a2.fm_line[np.isfinite(a2.fm_line)].mean(),
                        fm_nofirn_years=int((~np.isfinite(a2.fm_line)).sum()),
                        th_line_mean=a.th_line[np.isfinite(a.th_line) & a.th_line.notna()].mean(),
                        th_nofirn_years=int((a.th_line == np.inf).sum()),
                        ov_area=a2.ov_area.mean(),
                        gg_firn_fm_ice_km2=a2.gg_firn_fm_ice.mean(), fm_firn_gg_ice_km2=a2.fm_firn_gg_ice.mean()))
    return pd.DataFrame(res), ann


def zone(out, g, lo, hi):
    d = out[(out.band_mid >= lo) & (out.band_mid <= hi)].copy()
    d['decade'] = (d.year // 10) * 10
    d['w'] = d.rgi_area_km2
    rows = []
    for dec, a in d.groupby('decade'):
        a = a[a.mb_year_present.fillna(False).astype(bool)]
        if not len(a):
            continue
        w = a.w
        fm_frac = (a.fm_firn_frac * w).sum() / w[a.fm_firn_frac.notna()].sum()
        thk = a[a.fm_firn_frac >= 0.5]
        rows.append(dict(glacier=g, decade=f'{dec}s', zone=f'{lo}-{hi}',
                         fm_firn_frac=fm_frac,
                         fm_firn_m_firncells=a.fm_firn_m_firncells.mean(),
                         fm_firn_m_allcells=(a.fm_firn_mwe_allcells * 1000 / RHO_FIRN).mean(),
                         th_firn_frac=(a.th_firn * w).sum() / w[a.th_firn.notna()].sum() if a.th_firn.notna().any() else np.nan,
                         th_wet_depth_m=a.wet_depth_m.mean(), th_water_lost=a.water_lost_mwe.mean(),
                         z_perm_m=a.z_perm_m.mean(), mb_mwe=a.mb_mwe.mean()))
    return pd.DataFrame(rows), d


if __name__ == '__main__':
    od = f'{S}/band_tables_v2'; os.makedirs(od, exist_ok=True)
    allo, alls, allz, alla = [], [], [], []
    for g, gid in GL.items():
        out, top = band_tables(g, gid)
        out.to_csv(f'{od}/annual_bands_{g}_{gid}.csv', index=False, float_format='%.4f')
        allo.append(out)
        s, ann = summarise(out, g); alls.append(s); ann.insert(0, 'glacier', g); alla.append(ann)
        if g in ZONES:
            z, zd = zone(out, g, *ZONES[g]); allz.append(z)
            zd.to_csv(f'{od}/zone_annual_{g}.csv', index=False, float_format='%.4f')
    pd.concat(allo).to_csv(f'{od}/annual_bands_all.csv', index=False, float_format='%.4f')
    pd.concat(alla).to_csv(f'{od}/annual_firnline_all.csv', index=False, float_format='%.2f')
    S_ = pd.concat(alls); S_.to_csv(f'{od}/decade_summary.csv', index=False, float_format='%.3f')
    Z = pd.concat(allz); Z.to_csv(f'{od}/zone_decade.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
    print(S_.round(2).to_string(index=False))
    print(Z.round(2).to_string(index=False))
