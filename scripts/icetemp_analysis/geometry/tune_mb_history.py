"""Test 1c, part 2 (plan of 20 Oct): the mass-balance history correction for the thermal test. Per period (1942-1971 =
the thermal spin-up's reference period, 1972-1985, 1986-1999; none from 2000, the calibration period), one
precipitation factor matched to the observed winter balance and one temperature offset matched to the remaining annual
balance, on the 5 western GLAMOS glaciers observed since before 1972 (Aletsch, Allalin, Hohlaub, Schwarzberg; Gietro
for the annual balance only). Mass-balance runs of the 8 long series (GloGEM branch mb-history-diagnostic,
switch mb_hist_file), iterated; writes the file the thermal runs read.

Usage: python scripts/icetemp_analysis/geometry/tune_mb_history.py [N_ITER]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import past_mb_check as pm   # noqa: E402

G = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/geom')
TREE = '/home/jabeer/projects/glogemflow_development/GloGEM_wt_mbh'
BASE = G / 'mb_glamoslong_1941'
WEST = ('aletsch', 'allalin', 'hohlaub', 'schwarzberg', 'gietro')
PERIODS = ((1942, 1971), (1972, 1985), (1986, 1999))
SENS = -0.6        # m w.e. a-1 per K, first guess of the annual balance's temperature sensitivity


def write_file(path, P, dT, note):
    with open(path, 'w') as f:
        f.write(f'# mass-balance history correction, test 1c of 2026-10-20 ({note}): year  precipitation factor  dT [C]\n')
        for y in range(1941, 2026):
            k = next((i for i, (a, b) in enumerate(PERIODS) if a <= y <= b), 0 if y == 1941 else None)
            f.write(f'{y} {P[k] if k is not None else 1.0:.4f} {dT[k] if k is not None else 0.0:.4f}\n')


def run(tag, hist):
    d = G / 'mbh_tune' / tag
    (d / 'monthly/CentralEurope/calibration').mkdir(parents=True, exist_ok=True)
    for f in (BASE / 'monthly/CentralEurope/calibration').glob('*.dat'):
        (d / 'monthly/CentralEurope/calibration' / f.name).write_bytes(f.read_bytes())
    cfg = (BASE / 'config.pro').read_text().replace(f"dirres = '{BASE}/'", f"dirres = '{d}/'")
    (d / 'config.pro').write_text(cfg + f"mb_hist_file = '{hist}'\n")
    subprocess.run(f"cd {TREE} && GLOGEM_CONFIG={d}/config.pro bash -c \"echo '.r glogem' | nice -n 5 idl\" > {d}/idl.log 2>&1",
                   shell=True, check=True)
    return d / 'monthly/CentralEurope/PAST/PAST_original'


def stats(run_dir):
    pm.RUN = run_dir
    Ba, Bw = pm.model('Annual_Balance_sfc'), pm.model('Winter_balance_sfc')
    out = []
    for y0, y1 in PERIODS:
        bw_o, bw_m, ba_o, ba_m = [], [], [], []
        for name in WEST:
            o, gid = pm.glamos(name), pm.GLACIERS[name]
            yrs = [y for y in range(y0, y1 + 1) if y in o.index and y in Ba.columns]
            ba_o += list(o.Ba.loc[yrs]); ba_m += list(Ba.loc[gid, yrs])
            wy = [y for y in yrs if np.isfinite(o.Bw.get(y, np.nan))]
            bw_o += list(o.Bw.loc[wy]); bw_m += list(Bw.loc[gid, wy])
        out.append(dict(Bw_obs=np.mean(bw_o), Bw_mod=np.mean(bw_m), Ba_obs=np.mean(ba_o), Ba_mod=np.mean(ba_m)))
    return out


def main():
    n_iter = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    (G / 'mbh_tune').mkdir(exist_ok=True)
    P, dT = [1.0] * 3, [0.0] * 3
    s = stats(BASE / 'monthly/CentralEurope/PAST/PAST_original')
    sens = [SENS] * 3
    for it in range(n_iter):
        P = [P[k] * s[k]['Bw_obs'] / s[k]['Bw_mod'] for k in range(3)]
        hist = G / 'mbh_tune' / f'hist_it{it}a.dat'
        write_file(hist, P, dT, f'iteration {it}, precipitation step')
        sa = stats(run(f'it{it}a', hist))
        dT_new = [dT[k] + (sa[k]['Ba_mod'] - sa[k]['Ba_obs']) / -sens[k] for k in range(3)]
        hist = G / 'mbh_tune' / f'hist_it{it}b.dat'
        write_file(hist, P, dT_new, f'iteration {it}, temperature step')
        s = stats(run(f'it{it}b', hist))
        sens = [(s[k]['Ba_mod'] - sa[k]['Ba_mod']) / (dT_new[k] - dT[k]) if abs(dT_new[k] - dT[k]) > 0.05 else sens[k]
                for k in range(3)]
        dT = dT_new
        for k, (y0, y1) in enumerate(PERIODS):
            print(f'it {it} {y0}-{y1}: P {P[k]:.3f} dT {dT[k]:+.2f} | Bw mod {s[k]["Bw_mod"]:.2f} obs {s[k]["Bw_obs"]:.2f} | '
                  f'Ba mod {s[k]["Ba_mod"]:.2f} obs {s[k]["Ba_obs"]:.2f} | sens {sens[k]:.2f}', flush=True)
    write_file(G / 'mb_history_west.dat', P, dT, f'final after {n_iter} iterations')
    print('written', G / 'mb_history_west.dat')


if __name__ == '__main__':
    main()
