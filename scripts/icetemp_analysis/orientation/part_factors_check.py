"""Checks of part_factors.py before its tables are used (step 1c of the Step 5 plan).

  sun:     the vectorised June-September radiation (part_factors.Sun) equals orientation_eb.clear_sky_hourly at the
           same latitude and elevation (plane direct, flat direct, flat diffuse), to 1e-9 relative
  factor:  the energy-balance path of part_factors (psi0 anchor, eb_jjas) with each hole's own geometry and the band
           inputs of Stage 1 (band_forcing.csv) reproduces the Stage-1 factors (orient_tables/factors.csv, central)
           at every dT, to 1e-6 (the psi0 root is found to 0.02 W m-2 on sums in another order)
  inputs:  the logged band climate of the melth run (firnice_tairth/melth) against the Stage-1 rebuild
           (band_forcing.csv) at the borehole bands: air temperature and June-September degree-day melt, 1981-2020

Usage: PYTHONPATH=src python scripts/icetemp_analysis/orientation/part_factors_check.py MELTH_RUN_DIR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import orientation_eb as oe  # noqa: E402
import part_factors as pf  # noqa: E402
import write_orient_tables as wot  # noqa: E402

N_HOLES = 12


def check_sun(rng):
    worst = 0.0
    for _ in range(6):
        lat, z = rng.uniform(45.8, 46.6), rng.uniform(2500, 4600)
        slope, aspect = rng.uniform(0, 0.8, 5), rng.uniform(0, 2 * np.pi, 5)
        hor = rng.uniform(-0.05, 0.35, (5, 72))
        sun = pf.Sun(lat)
        beam, fdif, fdir = sun.beam(z)
        sdir = sun.direct(beam, slope, aspect, hor)
        for i in range(5):
            ref = oe.clear_sky_hourly(lat, z, slope[i], aspect[i], hor[i])
            for a, b in ((sdir[i], ref[0][pf.JJAS_DAYS]), (pf.hour_means(fdir), ref[1][pf.JJAS_DAYS]),
                         (pf.hour_means(fdif), ref[2][pf.JJAS_DAYS])):
                worst = max(worst, float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-9)))
    return worst


def check_factor(rng):
    w = oe.read_station(wot.INPUTS['jungfraujoch'])
    sig, dif, _ = oe.sky_weights(w)
    F = pd.read_csv(wot.INPUTS['band_forcing'])
    H = pd.read_csv(wot.INPUTS['holes'], dtype={'hole': str, 'rgi': str})
    comp = np.load(wot.INPUTS['horizons'])
    hor_of = {int(i): h for i, h in zip(comp['hole_id'], comp['hor'])}
    ref = pd.read_csv(pf.W3 / 'orient_tables/factors.csv', dtype={'hole': str, 'rgi': str})
    ref = ref[ref.strength == 'C'].set_index(['hole', 'dT']).F
    f = F[(F.year >= oe.Y0) & (F.year <= oe.Y1)]
    nyr = f.groupby('entity').year.nunique()
    dd_m = (f.groupby(['entity', 'month']).melt_dd.sum().unstack() * 1000.0).div(nyr, axis=0)
    sig_j = sig.reshape(oe.NY, 365)[:, pf.JJAS_DAYS, None]
    dif_j = dif.reshape(oe.NY, 365)[:, pf.JJAS_DAYS, None]
    worst, n = 0.0, 0
    ents = rng.choice(sorted(H.entity.unique()), N_HOLES, replace=False)
    for e in ents:
        s = H[H.entity == e]
        Th = oe.Forcing(oe.band_T(F, e), w, sig, dif).Th[:, pf.JJAS_DAYS, :]
        target = max(oe.jjas(dd_m.loc[e].reindex(range(1, 13)).fillna(0).to_numpy()), pf.DMIN)
        geo = []
        for _, h in s.iterrows():
            sd, fdr, fdf, _ = oe.clear_sky_hourly(h.lat, h.elevation, np.radians(h.slope_deg), np.radians(h.aspect_deg),
                                                  hor_of[int(h.hole)])
            S = sig_j * sd[pf.JJAS_DAYS][None] + dif_j * h.svf * fdf[pf.JJAS_DAYS][None]
            Sf = sig_j * fdr[pf.JJAS_DAYS][None] + dif_j * fdf[pf.JJAS_DAYS][None]
            geo.append((h.hole, S, Sf))
        p0 = brentq(lambda x: pf.eb_jjas(Th, geo[0][2], x) - target, -900, 400, xtol=0.02)
        for hole, S, Sf in geo:
            for dT in pf.DTS:
                Fv = pf.eb_jjas(Th + dT, S, p0) / max(pf.eb_jjas(Th + dT, Sf, p0), 1e-9)
                worst = max(worst, abs(Fv - ref[(hole, dT)]))
                n += 1
    return worst, n


def check_inputs(run):
    F = pd.read_csv(wot.INPUTS['band_forcing'], dtype={'rgi': str})
    rows = []
    for (e, gid, bm), d in F.groupby(['entity', 'rgi', 'band_model']):
        inp = pf.band_inputs(run, gid, [float(bm)])[float(bm)]
        Tm = oe.band_T(F, e)
        f = d[(d.year >= oe.Y0) & (d.year <= oe.Y1)]
        dd = f[f.month.between(6, 9)].melt_dd.sum() * 1000.0 / oe.NY
        rows.append(dict(entity=e, T_max_diff=float(np.max(np.abs(inp['Tm'] - Tm))), dd_logged=inp['dd_jjas'],
                         dd_rebuilt=dd))
    return pd.DataFrame(rows)


def main():
    run = Path(sys.argv[1])
    rng = np.random.default_rng(3)
    s = check_sun(rng)
    print(f'sun: largest relative difference to clear_sky_hourly {s:.1e}  {"PASS" if s < 1e-9 else "FAIL"}')
    fwor, n = check_factor(rng)
    print(f'factor: largest difference to the Stage-1 factors over {n} hole-dT values {fwor:.1e}  '
          f'{"PASS" if fwor < 1e-6 else "FAIL"}')
    I = check_inputs(run)
    I['dd_ratio'] = I.dd_logged / I.dd_rebuilt.where(I.dd_rebuilt > 1)
    print(f'inputs at {len(I)} borehole bands: air temperature, largest difference {I.T_max_diff.max():.3f} C; '
          f'June-September degree-day melt logged / rebuilt: median {I.dd_ratio.median():.3f}, '
          f'5-95 % {I.dd_ratio.quantile(.05):.3f}-{I.dd_ratio.quantile(.95):.3f}')
    print(I.sort_values('T_max_diff').tail(5).round(3).to_string(index=False))
    sys.exit(0 if s < 1e-9 and fwor < 1e-6 else 1)


if __name__ == '__main__':
    main()
