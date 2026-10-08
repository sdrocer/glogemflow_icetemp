"""Snow conduction test: GloGEM's thermal model on one ice band (snow_conduction_harness.pro)
against the exact steady answer.

The test snowpack: SNO m w.e. stacked from the surface, filling the gaps between nodes in turn,
each gap at the firn-curve density of the node above it (250, 300, 360 ... kg/m3). In steady
state the heat flow q is the same through every gap, so the temperature rises across a gap by
q x thickness / conductivity, with snow and ice in a gap adding their resistances.

Run from the repo root: python scripts/icetemp_analysis/snow_conduction_check.py <out_dir> [host]
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

GLOGEM = Path(__file__).resolve().parents[3] / 'GloGEM'
HARNESS = Path(__file__).with_name('snow_conduction_harness.pro')
CONFIG = (Path(__file__).resolve().parents[2]
          / 'data/bayescal/snowlog_w0/training_runs/config_bayescal_training.pro')
T_SURF, Q, K_ICE = -30.0, 0.3, 2.1
DENS_Z = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 19, 24, 29, 34, 39, 44, 49]          # settings.pro
DENS_REF = [250, 300, 360, 420, 480, 550, 580, 610, 640, 670, 700, 730, 760, 790, 820, 850, 880, 900]


def k_snow(rho):
    return min(0.024 - 1.23e-4 * rho + 2.5e-6 * rho ** 2, K_ICE)       # Calonne (2011), f_ins = 1


def exact(z, sno):
    """Exact steady temperature at the node depths z for the test snowpack."""
    T, rest = [T_SURF], sno
    for j in range(1, len(z)):
        dz = z[j] - z[j - 1]
        rho = np.interp(z[j - 1], DENS_Z, DENS_REF)
        s = min(rest / (rho / 1000), dz) if rest > 0 else 0.0
        rest -= s * rho / 1000
        T.append(T[-1] + Q * (s / k_snow(rho) + (dz - s) / K_ICE))
    return np.array(T)


def run(sno, implicit, out, host):
    env = (f"GLOGEM_CONFIG='{CONFIG}' HARNESS_SNO={sno} HARNESS_IMPLICIT={implicit} "
           f"HARNESS_OUT='{out}'")
    cmd = f"cd {GLOGEM} && {env} idl -e '.r {HARNESS}'"
    subprocess.run(['ssh', '-o', 'BatchMode=yes', host, cmd], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return np.loadtxt(out)


def main():
    out_dir = Path(sys.argv[1]); out_dir.mkdir(parents=True, exist_ok=True)
    host = sys.argv[2] if len(sys.argv) > 2 else 'vierzack03'
    worst = 0.0
    for implicit in ('y', 'n'):
        for sno in (0.0, 0.15, 0.30, 0.45):
            d = run(sno, implicit, out_dir / f'col_{implicit}_{sno}.txt', host)
            z, T = d[:, 0], d[:, 1]
            Te = exact(z, sno)
            err = np.abs(T - Te).max()
            worst = max(worst, err)
            k = int(np.searchsorted(z, 3.0))
            print(f'implicit={implicit} snow {sno:.2f} m w.e.: max |model - exact| {err:.5f} C '
                  f'(3 m: model {T[k]:.3f}, exact {Te[k]:.3f}; bed: {T[-1]:.3f} / {Te[-1]:.3f})')
    print(f'worst {worst:.5f} C')


if __name__ == '__main__':
    main()
