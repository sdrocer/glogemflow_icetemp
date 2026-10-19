"""Firn-model thickness against the thermal model's wet depth on its firn bands, per glacier and decade."""
import pandas as pd, numpy as np, os
S = os.environ.get('FIRNMODEL_DIR', '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel')
d = pd.read_csv(f'{S}/band_tables_v2/annual_bands_all.csv')
d = d[(d.glamos_cells >= 3) & (d.rgi_area_km2 > 0) & (d.th_firn == 1) & d.mb_year_present.fillna(False).astype(bool)]
d['fm_m_all'] = d.fm_firn_mwe_allcells * 1000 / 600.
d['decade'] = (d.year // 10 * 10).astype(str) + 's'
def agg(a):
    w = a.rgi_area_km2
    wa = lambda c: (a[c] * w)[a[c].notna()].sum() / w[a[c].notna()].sum()
    return pd.Series(dict(th_firn_area_km2=w.sum() / a.year.nunique(), fm_firn_frac=wa('fm_firn_frac'),
                          fm_thick_m=wa('fm_m_all'), wet_depth_m=wa('wet_depth_m'), z_perm_m=wa('z_perm_m')))
r = d.groupby(['glacier', 'decade']).apply(agg).reset_index()
r.to_csv(f'{S}/band_tables_v2/thermal_firn_bands_thickness_vs_wetdepth.csv', index=False, float_format='%.3f')
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)
print(r.round(2).to_string(index=False))
