"""Check of GloGEM's firnice_orient = 'table' (orient_harness.pro), run in IDL and judged here.

  read:   the dT line, and table lines to the band within 5 m of their elevation (also off the band centre); other
          bands F = 1; on the production grid -4..+6 K and on an uneven grid
  tref:   the 1981-2020 mean of each calendar month at the cell plus the lapse rate to the band
  factor: melt with the table / melt without = F interpolated linearly at the month's anomaly, held at the ends,
          to 1e-12; an all-ones table gives the same melt bit for bit
  state:  the mass-balance variables are untouched
  errors: a glacier without a table, and a band given twice, stop the run with a message

Run from the repo root: python scripts/icetemp_analysis/orientation/orient_check.py <out_dir> [glogem_tree] [host]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

HARNESS = Path(__file__).with_name('orient_harness.pro')
CONFIG = (Path(__file__).resolve().parents[3]
          / 'data/bayescal/ce25_pred_w0/training_runs/config_bayescal_training.pro')
GRIDS = dict(random=np.arange(-4.0, 7.0), uneven=np.array([-4.0, -2.5, -1.0, 0.0, 0.5, 2.0, 6.0]),
             ones=np.arange(-4.0, 7.0))
NB, ELEV = 30, 2600.0 + 70.0 * np.arange(30)


def run(tree, table_dir, out, host):
    cmd = (f"cd {tree} && GLOGEM_CONFIG='{CONFIG}' HARNESS_DIR='{table_dir}' HARNESS_OUT='{out}' "
           f"idl -e '.r {HARNESS}' > '{out}.log' 2>&1")
    subprocess.run(['ssh', '-o', 'BatchMode=yes', host, cmd], check=True)
    if not Path(out).exists():
        raise SystemExit(f'no output, see {out}.log')
    lines = Path(out).read_text().split('\n')
    i_f, i_m = lines.index('F'), lines.index('months')
    tref = np.loadtxt(lines[1:i_f])
    F = np.loadtxt(lines[i_f + 1:i_m])
    M = np.loadtxt(lines[i_m + 1:])
    return tref, F, M


def stops(tree, table_dir, out, host, text):
    cmd = (f"cd {tree} && GLOGEM_CONFIG='{CONFIG}' HARNESS_DIR='{table_dir}' HARNESS_OUT='{out}' "
           f"idl -e '.r {HARNESS}' > '{out}.log' 2>&1")
    subprocess.run(['ssh', '-o', 'BatchMode=yes', host, cmd], check=False)
    log = Path(f'{out}.log').read_text()
    return text in log and 'months' not in (Path(out).read_text() if Path(out).exists() else '')


def expected_tref():
    nt = 60 * 12
    ryear = 1961 + np.arange(nt) // 12
    rmon = np.arange(nt) % 12 + 1
    temp = -2 - 9 * np.cos(2 * np.pi * (rmon - 1) / 12) + 0.03 * (ryear - 1990) + 0.7 * np.sin(1.7 * np.arange(nt))
    dtdz = -0.0058 - 0.0004 * np.cos(2 * np.pi * np.arange(12) / 12)
    sel = (ryear >= 1981) & (ryear <= 2020)
    clim = np.array([temp[sel & (rmon == m)].mean() for m in range(1, 13)])
    return clim[None, :] + (ELEV - 2500.0)[:, None] * dtdz[None, :]


def main():
    out = Path(sys.argv[1]).resolve()
    tree = sys.argv[2] if len(sys.argv) > 2 else '/home/jabeer/projects/glogemflow_development/GloGEM_wt_orient'
    host = sys.argv[3] if len(sys.argv) > 3 else 'vierzack03'
    rng = np.random.default_rng(1)
    rows = {3: 0.0, 10: 3.0, 11: -4.5, 22: 5.0, 29: 0.0}          # band: offset of the line's elevation [m]
    ok = True
    for name, DTS in GRIDS.items():
        tab = {j: (np.ones(len(DTS)) if name == 'ones' else rng.uniform(0.3, 1.8, len(DTS))) for j in rows}
        d = out / name
        d.mkdir(parents=True, exist_ok=True)
        (d / 'orient_99999.dat').write_text('# test table\ndT ' + ' '.join(f'{x:.1f}' for x in DTS) + '\n\n' + ''.join(
            f'{ELEV[j] + rows[j]:.1f} ' + ' '.join(f'{x:.15f}' for x in tab[j]) + '\n' for j in rows))
        tref, F, M = run(tree, d, out / f'{name}.txt', host)
        F_exp = np.ones((NB, len(DTS)))
        for j in rows:
            F_exp[j] = np.round(tab[j], 15)
        read_ok = np.allclose(F, F_exp, rtol=0, atol=1e-14)
        tref_err = float(np.max(np.abs(tref - expected_tref())))
        m, j, tl, mel_n, mel_o, changed = M.T
        m, j = m.astype(int), j.astype(int)
        dT = np.clip(tl - tref[j, m - 1], DTS[0], DTS[-1])
        f_exp = np.array([np.interp(x, DTS, F[jj]) for x, jj in zip(dT, j)])
        pos = mel_n > 0
        ratio_err = float(np.max(np.abs(mel_o[pos] / mel_n[pos] - f_exp[pos])))
        zero_ok = bool(np.all(mel_o[~pos] == 0))
        state_ok = bool(np.all(changed == 0))
        spread = (float(dT[pos].min()), float(dT[pos].max()), int(pos.sum()))
        good = read_ok and tref_err < 1e-9 and ratio_err < 1e-12 and zero_ok and state_ok
        if name == 'ones':
            ident = bool(np.array_equal(mel_o, mel_n))
            good &= ident
            print(f'ones: melt with the table identical to without: {ident}')
        print(f'{name}: table read {read_ok}, tref error {tref_err:.1e} C, factor error {ratio_err:.1e} over '
              f'{spread[2]} melting band-months (dT {spread[0]:+.2f} to {spread[1]:+.2f} K), no melt kept zero '
              f'{zero_ok}, mass balance untouched {state_ok}  {"PASS" if good else "FAIL"}')
        ok &= good
    d = out / 'missing'
    d.mkdir(parents=True, exist_ok=True)
    good = stops(tree, d, out / 'missing.txt', host, 'no table')
    print(f'missing table stops the run: {good}  {"PASS" if good else "FAIL"}')
    ok &= good
    d = out / 'twice'
    d.mkdir(parents=True, exist_ok=True)
    (d / 'orient_99999.dat').write_text('dT -1 1\n2950 1.1 1.2\n2951 0.9 0.8\n')
    good = stops(tree, d, out / 'twice.txt', host, 'given twice')
    print(f'a band given twice stops the run: {good}  {"PASS" if good else "FAIL"}')
    ok &= good
    print('ALL PASS' if ok else 'SOME CHECKS FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
