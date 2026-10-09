"""(d) bookkeeping: the model's snow at each T3 row's band and month (firnice_snowbh, identical in
every alps_s2 run), the band's firn/ice state (firnice_firnline: lowest thermal firn band of the
previous hydrological year, the mask the thermal model uses), and the snow depth in the model's
grid (z_mass: the snow fills the top layers at the firn-curve densities 250, 300, 360, ... kg/m3,
as firnice_temperature_model.pro's kgap loop and snow logbook do). Saves snow_rows.pkl."""
import sys, io, re, tarfile, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
ta = C.ta
pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 200)
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
RUN = C.REPO / 'data/bayescal/alps_s2_w0/training_runs/monthly/CentralEurope/PAST/firnice_temperature'
GIDS = ['00773', '00777', '01112', '01225', '01230', '02103', '03116']

# model grid and firn-curve densities (settings.pro)
fit_layers, fit_dzstep = [10, 10, 23], [1., 5., 20.]
dz = np.concatenate([np.full(n, s) for n, s in zip(fit_layers, fit_dzstep)])
znode = np.concatenate([[0.0], np.cumsum(dz[1:])])
dens_z = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 19, 24, 29, 34, 39, 44, 49], float)
dens_ref = np.array([250, 300, 360, 420, 480, 550, 580, 610, 640, 670, 700, 730, 760, 790, 820, 850, 880, 900], float)
fit_dens = np.clip(np.interp(znode, dens_z, dens_ref), dens_ref[0], 900.0)


def z_mass(sno):
    """Snow depth [m] in the model grid for sno [m w.e.] (snow logbook loop, firnice_temperature_model.pro)."""
    rest, zm = float(sno), 0.0
    for j in range(1, len(dz)):
        if rest <= 0:
            break
        rho = fit_dens[j - 1]
        d = min(rest / (rho / 1000.0), dz[j])
        rest -= d * rho / 1000.0
        zm += d
    return zm


# snowbh from design 0's archive (identical in all runs: md5 checked for designs 0, 13, 31)
tgz = ta.archive_path('alps_s2', 0)
snow = {}
with tarfile.open(tgz) as tar:
    for gid in GIDS:
        txt = tar.extractfile(f'firnice_snowbh_{gid}.dat').read().decode()
        df = pd.read_csv(io.StringIO(txt), sep=r'\s+', comment='#', names=['profile', 'band', 'year', 'month', 'sno'])
        df['cal_year'] = df.year - (df.month >= 10)
        df['ym'] = df.cal_year * 100 + df.month
        for (p, b), g in df.groupby(['profile', 'band']):
            snow[(gid, int(p))] = (int(b), dict(zip(g.ym, g.sno)))
# firn line: lowest thermal firn band at the end of each hydrological year
fl = {}
for gid in GIDS:
    df = pd.read_csv(RUN / f'firnice_firnline_{gid}.dat', sep=r'\s+', comment='#', names=['year', 'ela', 'aar', 'fl', 'nfi', 'nif'])
    fl[gid] = dict(zip(df.year, df.fl))


def is_ice(gid, band, ym):
    y, mth = ym // 100, ym % 100
    hy = y + (mth >= 10)
    thr = fl[gid][max(hy - 1, 1972)]
    return band < thr


bu = mod['band_used']
prof = {}
for e, g in enumerate(obs['glaciers']):
    f = bu.file[e]
    n = int(re.match(r'temp_ID(\d+)_', f).group(1))
    band, ser = snow[(g.glacier_id, n)]
    if abs(band - bu.band_used[e]) > 0.5:
        raise RuntimeError(f'{g.glacier_name}: snowbh band {band} vs profile band {bu.band_used[e]}')
    prof[e] = (g.glacier_id, band, ser)

T3 = obs['T3'].copy()
pm = ta.prev_month
rec = []
for r in T3.itertuples():
    gid, band, ser = prof[r.e]
    ym, yp = int(r.ym), pm(int(r.ym))
    s_m, s_p = ser.get(ym, np.nan), ser.get(yp, np.nan)
    ice_m, ice_p = is_ice(gid, band, ym), is_ice(gid, band, yp)
    rec.append(dict(sno_m=s_m, sno_p=s_p, ice_m=ice_m, ice_p=ice_p,
                    h_m=z_mass(s_m) if (ice_m and s_m > 0) else 0.0,
                    h_p=z_mass(s_p) if (ice_p and np.isfinite(s_p) and s_p > 0) else 0.0))
S = pd.concat([T3[['t3', 't2', 'e', 'entity', 'depth', 'ym', 'frac', 'count', 'firn_label', 'range_only', 'n_range']].reset_index(drop=True),
               pd.DataFrame(rec)], axis=1)
S['grp'] = ta.group_of(obs['T3'])
S['common'] = mod['common3']
S['snow_ice'] = S.ice_m & (S.sno_m > 0)                  # thermal step of the reading's month had snow on an ice band
S['snow_firn'] = ~S.ice_m & (S.sno_m > 0)
S['shallow'] = S.depth <= 3.0
# validate the reconstructed firn mask against snowclim's firn_frac (mean over the run, per band)
chk = []
with tarfile.open(tgz) as tar:
    for gid in GIDS:
        txt = tar.extractfile(f'firnice_snowclim_{gid}.dat').read().decode()
        sc = pd.read_csv(io.StringIO(txt), sep=r'\s+', comment='#',
                         names=['elev', 'month', 'firn_frac', 'sno', 'z_snow', 'z_mass', 'Ts', 'Ti', 'T1', 'T2', 'Ru', 'Rs'])
        for (g2, band, ser) in {prof[e] for e in prof if prof[e][0] == gid and True} if False else []:
            pass
        bands = sorted({prof[e][1] for e in prof if prof[e][0] == gid})
        for b in bands:
            ff = sc[sc.elev == b].firn_frac.mean()
            yms = [y * 100 + m for y in range(1971, 2026) for m in range(1, 13) if 197110 <= y * 100 + m <= 202509]
            rec_ff = np.mean([not is_ice(gid, b, ym) for ym in yms])
            chk.append((gid, b, round(float(ff), 3), round(float(rec_ff), 3)))
chk = pd.DataFrame(chk, columns=['gid', 'band', 'firn_frac_snowclim', 'firn_frac_reconstructed'])
print('firn mask check (fraction of months firn):'); print(chk.to_string(index=False))
print('largest difference:', float((chk.firn_frac_snowclim - chk.firn_frac_reconstructed).abs().max()))
# z_mass check against snowclim's monthly means for HEF 2615 m
c = S[S.common]
print('\nT3 common rows:', len(c), '; on model ice bands:', int(c.ice_m.sum()), '; with snow in the month on ice bands:', int(c.snow_ice.sum()),
      '; on firn bands with snow:', int(c.snow_firn.sum()))
sh = c[c.shallow]
print('T3 rows <= 3 m:', len(sh), ' of them on ice bands with snow in the month:', int(sh.snow_ice.sum()),
      ' (prev month only:', int((~sh.snow_ice & sh.ice_p & (sh.sno_p > 0)).sum()), '); on firn bands with snow:', int(sh.snow_firn.sum()))
print(sh.groupby(['grp', 'firn_label']).agg(rows=('t3', 'size'), ice_band=('ice_m', 'sum'), snow_ice=('snow_ice', 'sum'),
                                            snow_firn=('snow_firn', 'sum')).to_string())
print('\nshallow rows on ice bands with snow:')
print(sh[sh.snow_ice][['entity', 'depth', 'ym', 'frac', 'count', 'sno_p', 'sno_m', 'h_p', 'h_m', 'range_only']].round(3).to_string())
# T2 points <= 3 m with any snow month
t2c = obs['T2'][mod['common2']]
pts = sh.groupby('t2').agg(any_snow=('snow_ice', 'any'), cnt_snow=('count', lambda v: 0), n=('t3', 'size'))
pts['w_snow'] = sh[sh.snow_ice].groupby('t2')['count'].sum().reindex(pts.index).fillna(0)
pts['w_all'] = sh.groupby('t2')['count'].sum()
print('\nT2 points <= 3 m:', len(pts), ' with snow on an ice band in some reading month:', int(pts.any_snow.sum()),
      '; their readings under snow:', float(pts.w_snow.sum()), 'of', float(pts.w_all[pts.any_snow].sum()))
# deeper rows on snow-covered ice bands (all depths)
print('all T3 common rows on ice bands with snow in the month, by depth band:',
      c[c.snow_ice].groupby(ta.depth_band(c[c.snow_ice].depth)).size().to_dict())
pickle.dump(dict(S=S, z_mass_table={s: z_mass(s) for s in (0.1, 0.25, 0.5, 0.7, 1.0)}, firn_check=chk), open(C.OUT / 'snow_rows.pkl', 'wb'))
print('z_mass(0.1, 0.25, 0.5, 0.7, 1.0 m w.e.) =', {s: round(z_mass(s), 3) for s in (0.1, 0.25, 0.5, 0.7, 1.0)})
