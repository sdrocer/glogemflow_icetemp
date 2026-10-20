"""Step 1a of the plan of 20 Oct (development_history/plans/2026-10-20_geometry_and_step5_plan.md): does thin model ice
line up with the warm misfits? Per borehole of the scored set (CentralEurope_fair without 01161, one point per borehole
and year):
- the model's band ice thickness (2003 inventory, fixed in the thermal runs; firnice_bed_<id>.dat of the stage-2 run);
- the drilled depth and whether the hole reached the bed (glenglat);
- the band thickness in the hole's measurement year: the 2003 thickness plus the ice lost between that year and 2003,
  the glacier's volume change (mass balance of the stage-2 settings, week3/geom/mb_1972_2025) spread over the bands
  with the Δh curve of glacier_retreat.pro (Huss et al. 2010), as GloGEM's Δh parameterisation would;
- the thinning 1972-2003 of the band (the dh of test 1b);
- the misfit of FRHS (each massif at its chosen f_ref) at the hole's own points below 15 m, and at its deepest points.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/geometry/thickness_vs_misfit.py [OUT_DIR]
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import forcing_fair_test as ft
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
MB = W3 / 'geom/mb_1972_2025/monthly/CentralEurope/PAST/PAST_original'
BANDS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope')
VARIANT, INVENTORY, Y0 = 'FRHS', 2003, 1972
DEEP, DEEPEST = 15.0, 10.0          # points below 15 m; 'deepest' = within 10 m of the hole's deepest point
DH_SIZE = (5.0, 20.0)               # km2, settings.pro dh_size
FREF_UNSCORED = 0.85                # massifs without scored points have no chosen f_ref


def table(var):
    d = pd.read_csv(MB / f'centraleurope_{var}_r1_CentralEurope_fair.dat', sep=r'\s+', dtype={'ID': str}).set_index('ID')
    d.columns = d.columns.astype(int)
    return d


def bands(gid):
    b = np.atleast_2d(np.loadtxt(BANDS / f'{gid}.dat', skiprows=5))
    return pd.DataFrame(dict(elev=((b[:, 1] + b[:, 2]) / 2).round().astype(int), area=b[:, 3], thick=b[:, 4]))


def dh_shape(bd):
    """glacier_retreat.pro: the normalised thickness change per band (1 at the terminus), 0 where there is no ice."""
    ice = bd.thick > 0
    e = bd.elev[ice].to_numpy()
    hr = (e.max() - e) / (e.max() - e.min())
    ta = bd.area.sum()
    if ta > DH_SIZE[1]:
        dh = (hr - 0.02) ** 6 + 0.12 * (hr - 0.02)
    elif ta > DH_SIZE[0]:
        dh = (hr - 0.05) ** 4 + 0.19 * (hr - 0.05) + 0.01
    else:
        dh = (hr - 0.30) ** 2 + 0.60 * (hr - 0.30) + 0.09
    out = np.zeros(len(bd))
    out[ice.to_numpy()] = np.clip(dh, 0, 1)
    return out


def ice_lost(gid, y_from, y_to, B, A, bd, shape):
    """Ice thickness lost per band [m] between the end of hydrological year y_from and the end of y_to (negative:
    gained), the glacier-wide balance spread with the Δh curve on the 2003 geometry."""
    yrs = [y for y in B.columns if y_from < y <= y_to]
    dvol = float((B.loc[gid, yrs] * A.loc[gid, yrs]).sum()) / 0.9 * 1e6     # m3 of ice, negative = loss
    return -dvol / float((shape * bd.area * 1e6).sum()) * shape


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'geom'
    out.mkdir(parents=True, exist_ok=True)
    res = json.load(open(W3 / 'stage2/stage2_hole_year.json'))
    chosen = res['summary'][VARIANT]['chosen_fref']
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    bh = ts.load_borehole_table().set_index('id')
    B, A = table('Annual_Balance_sfc'), table('Area')
    outputs = {}
    rows = []
    for g in dh.calibration_glaciers:
        if not g.glacier_id or g.glacier_id in ft.SEALED or g.glacier_id not in ft.MASSIF_OF:
            continue
        fr = chosen.get(ft.MASSIF_OF[g.glacier_id], FREF_UNSCORED)
        run = ft.run_name('s2', VARIANT, fr)
        if run not in outputs:
            outputs[run] = pickle.load(open(ft.output_path(run, 'hole_year'), 'rb'))
        o = outputs[run].get(g.glacier_name)
        if o is None:
            continue
        mod = np.asarray(o[1], float)[:len(g.depths)]
        bed = pd.read_csv(ft.BAYESCAL / run / ft.FIRNICE / f'firnice_bed_{g.glacier_id}.dat', sep=r'\s+', comment='#',
                          header=None, names=['elev', 'year', 'bed', 'T', 'Tmin', 'Tmax', 'pmp', 'bottom', 'Tb'])
        bed_band = bed.loc[(bed.elev - g.elevation_band).abs().idxmin()]
        bd = bands(g.glacier_id)
        shape = dh_shape(bd)
        k = int((bd.elev - g.elevation_band).abs().idxmin())
        for hole in sorted(set(g.point_hole)):
            sel = np.flatnonzero(g.point_hole == hole)
            d, ob, md, yr = g.depths[sel], g.T_obs[sel], mod[sel], g.point_year[sel]
            year = int(np.median(yr))
            lost_y = ice_lost(g.glacier_id, year, INVENTORY, B, A, bd, shape)[k] if year <= INVENTORY else \
                -ice_lost(g.glacier_id, INVENTORY, year, B, A, bd, shape)[k]
            deep = (d > DEEP) & np.isfinite(md)
            deepest = deep & (d >= d.max() - DEEPEST)
            h = bh.loc[str(hole)] if str(hole) in bh.index else None
            rows.append(dict(
                hole=hole, site=g.glacier_name.split('@')[0], entity=g.glacier_name, rgi=g.glacier_id,
                massif=ft.MASSIF_OF[g.glacier_id], band=g.elevation_band, year=year, fref=fr,
                drilled=float(h['depth']) if h is not None and pd.notna(h['depth']) else np.nan,
                to_bed=bool(str(h['to_bed']).lower() == 'true') if h is not None else False,
                deepest_reading=float(d.max()), model_thick_2003=float(bed_band.bed),
                model_column=float(bed_band.bottom), band_thick_inventory=float(bd.thick[k]),
                thick_meas_year=float(bd.thick[k] + lost_y), dh_1972_2003=float(ice_lost(g.glacier_id, Y0 - 1, INVENTORY,
                                                                                         B, A, bd, shape)[k]),
                n_deep=int(deep.sum()),
                bias_deep=float(np.mean(md[deep] - ob[deep])) if deep.any() else np.nan,
                rmse_deep=float(np.sqrt(np.mean((md[deep] - ob[deep]) ** 2))) if deep.any() else np.nan,
                bias_deepest=float(np.mean(md[deepest] - ob[deepest])) if deepest.any() else np.nan,
                obs_deepest=float(np.mean(ob[deepest])) if deepest.any() else np.nan))
    T = pd.DataFrame(rows)
    T['ratio_model_drilled'] = T.model_thick_2003 / T.drilled
    T['ratio_measyear_drilled'] = T.thick_meas_year / T.drilled
    T.to_csv(out / 'thickness_vs_misfit.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 250, 'display.max_rows', 300)
    cols = ['hole', 'entity', 'year', 'drilled', 'to_bed', 'deepest_reading', 'model_thick_2003', 'thick_meas_year',
            'dh_1972_2003', 'n_deep', 'bias_deep', 'bias_deepest']
    print(T[T.n_deep > 0][cols].round(2).to_string(index=False))
    tb = T[T.to_bed & (T.n_deep > 0)]
    print(f'\nto-bed holes with deep points: {len(tb)}; median model/drilled {tb.ratio_model_drilled.median():.2f}, '
          f'measurement-year/drilled {tb.ratio_measyear_drilled.median():.2f}')
    ok = T[(T.n_deep > 0) & np.isfinite(T.drilled)]
    print('rank correlation, deep bias against model/drilled thickness:',
          round(ok.bias_deep.corr(ok.ratio_model_drilled, method='spearman'), 2), f'(n={len(ok)})')


if __name__ == '__main__':
    main()
