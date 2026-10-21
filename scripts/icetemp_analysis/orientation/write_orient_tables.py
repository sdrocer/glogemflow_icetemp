"""Orientation factor tables for GloGEM (firnice_orient = 'table'): Stage 1 of the orientation test (Step 5 plan,
revised 21 Oct 2026, step 1b), F at the borehole bands for air temperature anomalies dT = -3..+3 K, three strengths:

  C  central:  c 10 W m-2 K-1, Jungfraujoch's daily direct beam spread over the day
  W  weaker:   c 20 (more of the melt energy follows the air temperature)
  S  stronger: c 10, the day's direct beam in sunny spells at full clear-sky strength nearest noon
  one          every F = 1 (the gate: must reproduce FRHS byte for byte)

Albedo 0.80, T_tip 1 C (orientation_eb.py). psi0 per band and strength: June-September melt of flat open ground equals
the band's degree-day thermal melt (FRHS, 1981-2020, at least 30 mm); it is kept when the hourly air temperature is
shifted by dT. F = June-September melt of the hole's plane over that of flat open ground, mean over the holes of the
outline's band (both sites where two share one). Bands without holes keep F = 1 in GloGEM.

Inputs: the orientation investigation (week3/orientation_investigation/): evidence/band_forcing.csv,
evidence/hole_sun.csv, data/calibration_holes_components.npz (horizons), MeteoSwiss Jungfraujoch daily data.
Check: F reproduces the investigation's designer/revise_F.csv, revise_F_dT.csv and revise_sens.csv.
Writes OUT/<C|W|S|one>/orient_<id>.dat, OUT/factors.csv (per hole) and OUT/tables.json; tables that exist with
other numbers are not overwritten.

Usage: PYTHONPATH=src python scripts/icetemp_analysis/orientation/write_orient_tables.py [OUT]
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_eb as oe  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
ORI = W3 / 'orientation_investigation'
INPUTS = dict(band_forcing=ORI / 'evidence/band_forcing.csv', holes=ORI / 'evidence/hole_sun.csv',
              horizons=ORI / 'data/calibration_holes_components.npz',
              jungfraujoch=ORI / 'downloads/meteoswiss_daily/ogd-smn_jun_d_historical.csv')
DTS = (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0)
STRENGTHS = {'C': dict(c=10.0, spells=False), 'W': dict(c=20.0, spells=False), 'S': dict(c=10.0, spells=True)}
ALB, TTIP, DMIN = 0.80, 1.0, 30.0


def band_of(entity):
    return int(entity.split('@')[1].split('~')[0])


def factors():
    w = oe.read_station(INPUTS['jungfraujoch'])
    sig, dif, fd_st = oe.sky_weights(w)
    spells = oe.concentrate(sig, fd_st)
    F = pd.read_csv(INPUTS['band_forcing'])
    H = pd.read_csv(INPUTS['holes'], dtype={'hole': str, 'rgi': str})
    comp = np.load(INPUTS['horizons'])
    hor_of = {int(i): h for i, h in zip(comp['hole_id'], comp['hor'])}
    f = F[(F.year >= oe.Y0) & (F.year <= oe.Y1)]
    nyr = f.groupby('entity').year.nunique()
    dd_m = (f.groupby(['entity', 'month']).melt_dd.sum().unstack() * 1000.0).div(nyr, axis=0)
    rows = []
    for e, s in H.groupby('entity'):
        Tm = oe.band_T(F, e)
        dd = dd_m.loc[e].reindex(range(1, 13)).fillna(0).to_numpy()
        target = max(oe.jjas(dd), DMIN)
        for name, o in STRENGTHS.items():
            fo = oe.Forcing(Tm, w, sig, dif, spells if o['spells'] else None)
            geo = []
            for _, h in s.iterrows():
                sd, fdr, fdf, _ = oe.clear_sky_hourly(h.lat, h.elevation, np.radians(h.slope_deg),
                                                      np.radians(h.aspect_deg), hor_of[int(h.hole)])
                geo.append((h.hole, fo.shortwave(sd, fdf, h.svf), fo.shortwave(fdr, fdf, 1.0)))
            p0 = brentq(lambda x: oe.jjas(oe.monthly_melt(fo.Th, geo[0][2], x, o['c'], ALB, TTIP)) - target,
                        -900, 400, xtol=0.02)
            for dT in DTS:
                Th = fo.Th + dT
                for hole, Ss, Sf in geo:
                    m_s = oe.jjas(oe.monthly_melt(Th, Ss, p0, o['c'], ALB, TTIP))
                    m_f = oe.jjas(oe.monthly_melt(Th, Sf, p0, o['c'], ALB, TTIP))
                    rows.append(dict(entity=e, rgi=s.rgi.iloc[0], band=band_of(e), hole=hole, strength=name, dT=dT,
                                     psi0=p0, dd_jjas=oe.jjas(dd), melt_slope=m_s, melt_flat=m_f,
                                     F=m_s / max(m_f, 1e-9)))
        print(e, {n: round(float(np.mean([r['F'] for r in rows if r['entity'] == e and r['strength'] == n
                                           and r['dT'] == 0.0])), 3) for n in STRENGTHS}, flush=True)
    return pd.DataFrame(rows)


def check(R):
    """Largest difference to the investigation's own factors (same construction)."""
    E = R.groupby(['entity', 'strength', 'dT']).F.mean()
    d = ORI / 'designer'
    a = pd.read_csv(d / 'revise_F.csv').set_index('entity')
    b = pd.read_csv(d / 'revise_F_dT.csv')
    c = pd.read_csv(d / 'revise_sens.csv').set_index('entity')
    diffs = {'C, dT 0': max(abs(E[(e, 'C', 0.0)] - a.at[e, 'F_c10']) for e in a.index),
             'W, dT 0': max(abs(E[(e, 'W', 0.0)] - a.at[e, 'F_c20']) for e in a.index),
             'S, dT 0': max(abs(E[(e, 'S', 0.0)] - c.at[e, 'F_JFJconc']) for e in c.index)}
    for cc, name in ((10.0, 'C'), (20.0, 'W')):
        bb = b[(b.c == cc) & (b.dT != 0.0)]
        diffs[f'{name}, dT +-1, +-2'] = max(abs(E[(r.entity, name, r.dT)] - r.F) for r in bb.itertuples())
    return {k: float(v) for k, v in diffs.items()}


def table_text(rgi, B, name):
    head = (f'# orientation factor F of the thermal melt, glacier {rgi}, strength {name} '
            f'(write_orient_tables.py, Step 5 plan of 21 Oct 2026, step 1b)\n'
            f'# band elevation [m]  F at dT = ' + ', '.join(f'{d:+.0f}' for d in DTS) + ' K\n')
    lines = [f'{band:d} ' + ' '.join(f'{x:.6f}' for x in vals) + '\n' for band, vals in B]
    return head + ''.join(lines)


def write(path, text):
    if path.exists() and path.read_text() != text:
        raise SystemExit(f'{path} exists with other numbers; use another folder')
    path.write_text(text)


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'orient_tables'
    R = factors()
    diffs = check(R)
    print('largest difference to the investigation:', diffs, flush=True)
    if max(diffs.values()) > 1e-6:
        raise SystemExit('the factors do not reproduce the investigation')
    band_F = R.groupby(['rgi', 'band', 'strength', 'dT']).F.mean()
    for name in list(STRENGTHS) + ['one']:
        (out / name).mkdir(parents=True, exist_ok=True)
        for rgi in sorted(R.rgi.unique()):
            bands = sorted(R[R.rgi == rgi].band.unique())
            src = name if name != 'one' else 'C'
            B = [(b, [1.0] * len(DTS) if name == 'one' else [band_F[(rgi, b, src, d)] for d in DTS]) for b in bands]
            write(out / name / f'orient_{rgi}.dat', table_text(rgi, B, name))
    R.to_csv(out / 'factors.csv', index=False, float_format='%.6f')
    at0 = R[R.dT == 0.0].groupby(['rgi', 'band', 'strength']).F.mean().unstack()
    meta = dict(constants=dict(albedo=ALB, T_tip=TTIP, min_flat_melt_mm=DMIN, dT=DTS,
                               strengths={k: dict(v) for k, v in STRENGTHS.items()}, weather='Jungfraujoch 1981-2020',
                               window='June-September', reference='flat open ground'),
                inputs={k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                        for k, p in INPUTS.items()},
                n_bands=int(len(at0)), n_holes=int(R.hole.nunique()), n_glaciers=int(R.rgi.nunique()),
                check_against_investigation=diffs,
                F_at_dT0={k: dict(min=float(at0[k].min()), median=float(at0[k].median()), max=float(at0[k].max()))
                          for k in STRENGTHS})
    write(out / 'tables.json', json.dumps(meta, indent=1) + '\n')
    print(json.dumps({k: v for k, v in meta.items() if k != 'inputs'}, indent=1))


if __name__ == '__main__':
    main()
