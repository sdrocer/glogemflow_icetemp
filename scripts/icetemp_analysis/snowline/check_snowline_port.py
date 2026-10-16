"""Gate 2 of the snowline-calibration port (plan of 2026-10-16): Python replicas of the IDL pieces.

  units   compare_snowline_rmse and snowline_p10 (GloGEM functions/) on synthetic data, run through IDL
  spline  her minimum selection (natural cubic spline through the scanned RMSE, 1000-point grid),
          fed with the RMSE curves of a scan log (*_snowline.dat), against the c_prec it chose

Usage: python check_snowline_port.py units GLOGEM_DIR
       python check_snowline_port.py spline LOG.dat [LOG2.dat ...]
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))


def rmse_py(gid, ids, year, doy, sla, sig, model, year0, start_doy=274):
    """compare_snowline_rmse: day 0 of `model` is day start_doy of year0-1 (365-day calendar)."""
    m = ids == gid
    k = 365 * (year[m].astype(int) - year0 + 1) + doy[m].astype(int) - start_doy
    ok = (k >= 0) & (k < len(model))
    if not ok.any():
        return np.nan, 0, np.nan
    e2 = (sla[m][ok] - model[k[ok]]) ** 2
    w = 1.0 / sig[m][ok]
    return np.sqrt(np.sum(e2 * w) / np.sum(w)), int(ok.sum()), np.sqrt(e2.mean())


def p10_py(sno, gl, area, elev, step, noval=-9999.0):
    ig = (gl != noval) & (area > 0)
    if not ig.any():
        return noval
    js = np.where((sno > 0) & ig)[0]
    if len(js) == 0:
        return elev[ig].max() + step / 2
    ca = np.cumsum(area[js]) / area[js].sum()
    k = int(np.argmax(ca >= 0.1))
    ca0 = ca[k - 1] if k > 0 else 0.0
    return elev[js[k]] - step / 2 + step * (0.1 - ca0) / (ca[k] - ca0)


def spline_min(x, y):
    """Her selection: natural cubic spline (IDL SPL_INIT default), minimum on findgen(1000)/999 grid."""
    if len(x) == 1:
        return x[0]
    xf = (np.arange(1000, dtype=np.float32) / np.float32(999.0)).astype(float) * (x[-1] - x[0]) + x[0]
    yf = CubicSpline(x, y, bc_type='natural')(xf) if len(x) > 2 else np.interp(xf, x, y)
    return xf[int(np.argmin(yf))]


def units(glogem_dir):
    rng = np.random.default_rng(16102026)
    tmp = tempfile.mkdtemp(prefix='snowport_')
    ny, year0 = 20, 2000
    model = 2500 + 600 * np.sin(np.arange(365 * ny) / 58.0) + rng.normal(0, 30, 365 * ny)
    n = 400
    ids = rng.choice(['01112', '02216', '00773'], n)
    year = rng.integers(1998, 2022, n).astype(float)
    doy = rng.integers(1, 366, n).astype(float)
    sla = rng.normal(2700, 300, n)
    sig = rng.uniform(80, 160, n)
    np.savetxt(f'{tmp}/model.txt', model, fmt='%.6f')
    with open(f'{tmp}/obs.dat', 'w') as f:
        f.write('header\n')
        for i in range(n):
            f.write(f'RGI70-11.{ids[i]} {sla[i]:.6f} {sig[i]:.6f} {int(year[i])} {int(doy[i])}\n')
    cases = [(g, s) for g in ('01112', '02216', '99999') for s in (274, 121)]
    # bands for p10: 40 bands of 10 m, some ice-free, snow on varying sets
    nb, step = 40, 10.0
    elev = 2000 + 5 + step * np.arange(nb)
    p10_cases = []
    for c in range(6):
        area = rng.uniform(0, 0.3, nb)
        gl = np.where(rng.random(nb) < 0.9, elev, -9999.0)
        area[gl == -9999.0] = 0
        sno = np.where(rng.random(nb) < [0.0, 1.0, 0.5, 0.2, 0.8, 0.5][c], rng.uniform(0.01, 2, nb), 0.0)
        p10_cases.append((sno, gl, area))
    with open(f'{tmp}/bands.txt', 'w') as f:
        for sno, gl, area in p10_cases:
            for arr in (sno, gl, area):
                f.write(' '.join(f'{v:.8f}' for v in arr) + '\n')
    cg = ','.join(f"'{g}'" for g, _ in cases)
    cs = ','.join(str(s) for _, s in cases)
    el = ','.join(f'{v:.1f}d' for v in elev)
    pro = f"""compile_opt idl2
!path = '{glogem_dir}/functions:{glogem_dir}/procedures/read:' + !path
read_snowlinedata, '{tmp}/obs.dat', gid, yo, dy, so, ss
model = dblarr({365 * ny}) & openr, 1, '{tmp}/model.txt' & readf, 1, model & close, 1
openw, 2, '{tmp}/idl_out.txt'
cg = [{cg}] & cs = [{cs}]
for i = 0, n_elements(cg) - 1 do begin
  r = compare_snowline_rmse(gid, cg[i], yo, dy, so, ss, model, {year0}, start_doy=cs[i], n_used=nu, rmse_unweighted=ru)
  printf, 2, cg[i], cs[i], r, nu, ru, format='(a,i5,e25.16,i8,e25.16)'
endfor
elev = [{el}]
b = dblarr({nb}) & openr, 1, '{tmp}/bands.txt'
for c = 0, {len(p10_cases) - 1} do begin
  readf, 1, b & sno = b & readf, 1, b & gl = b & readf, 1, b & area = b
  printf, 2, 'p10', snowline_p10(sno, gl, area, elev, {step}, -9999.d), format='(a,e25.16)'
endfor
close, 1 & close, 2
end
"""
    with open(f'{tmp}/run.pro', 'w') as f:
        f.write(pro)
    subprocess.run(f"cd {tmp} && echo '.run run.pro' | idl -quiet", shell=True, check=True, capture_output=True)
    out = open(f'{tmp}/idl_out.txt').read().split('\n')
    ok = True
    for line, (g, s) in zip(out[:len(cases)], cases):
        p = line.split()
        r_idl, nu_idl = float(p[2]), int(p[3])
        r_py, nu_py, _ = rmse_py(g, ids, year, doy, sla, sig, model, year0, s)
        same = (nu_idl == nu_py) and ((np.isnan(r_idl) and np.isnan(r_py)) or abs(r_idl - r_py) < 1e-8 * max(1, r_py))
        ok &= same
        print(f'rmse {g} start {s}: IDL {r_idl:.10f} ({nu_idl} obs)  Python {r_py:.10f} ({nu_py})  {"ok" if same else "DIFFERS"}')
    for line, (sno, gl, area) in zip(out[len(cases):len(cases) + len(p10_cases)], p10_cases):
        v_idl, v_py = float(line.split()[1]), p10_py(sno, gl, area, elev, step)
        same = abs(v_idl - v_py) < 1e-6
        ok &= same
        print(f'p10: IDL {v_idl:.6f}  Python {v_py:.6f}  {"ok" if same else "DIFFERS"}')
    print('units:', 'PASS' if ok else 'FAIL')
    return ok


def spline(logs):
    ok = True
    for fn in logs:
        t = pd.read_csv(fn, sep=r'\s+', dtype={'ID': str})
        cols = [c for c in t.columns[7:]]
        xs = np.array([float(c) for c in cols])
        for _, r in t[t.Status == 2].iterrows():
            n = int(r.Nscan)
            x, y = xs[:n], r[cols].values[:n].astype(float)
            c_py = spline_min(x, y)
            same = abs(c_py - r.CprecMin) <= 1.5 * (x[-1] - x[0]) / 999 + 1e-4   # one step of the 1000-point grid
            ok &= same
            print(f'{os.path.basename(fn)} {r.ID}: log {r.CprecMin:.3f}  Python {c_py:.3f}  {"ok" if same else "DIFFERS"}')
    print('spline:', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    mode = sys.argv[1]
    res = units(sys.argv[2]) if mode == 'units' else spline(sys.argv[2:])
    sys.exit(0 if res else 1)
