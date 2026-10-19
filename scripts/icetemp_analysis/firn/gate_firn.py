"""Gate of firnice_perm_limit = 'firn' (plan 2026-10-14_firn_ledger_plan.md): 1 default byte-identical, 2 mass-balance
files identical, 3 IDL ledger = Python reference (driven by the logged band balance), 4 wet depth <= ledger depth (+ this
year's snow), 5 the zone bands against the firn model (reported)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import firn_ledger as fl
B = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/bayescal')
SUB = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature'
REF = 'fair_FRH_r0.85_i1.017_a1.00_w0'
GATE, NEW = sys.argv[1], sys.argv[2]   # e.g. firngate_FRH_r0.85_i1.017_a1.00_w0 firn_FRHL_r0.85_i1.017_a1.00_w0
ok = {}
fa, fg, fn = ({p.name: p for p in (B / r / SUB).glob('*.dat')} for r in (REF, GATE, NEW))
d = [n for n in sorted(set(fa) & set(fg)) if fa[n].read_bytes() != fg[n].read_bytes()]
ok['1 default byte-identical'] = not d and set(fa) == set(fg)
print(f'1: {len(fa)} vs {len(fg)} files, differ {d[:4]}, only ref {sorted(set(fa) - set(fg))[:3]}, only gate {sorted(set(fg) - set(fa))[:3]}')
mb = [n for n in fa if n.startswith(('firnice_snowmay', 'firnice_snowbh', 'firnice_firnline'))]
d = [n for n in mb if fa[n].read_bytes() != fn[n].read_bytes()]
ok['2 mass-balance files identical'] = not d
print(f'2: {len(mb)} MB-derived files, differ {d[:4]}')
# 3: ledger against the Python reference, driven by the logged balance (column b if present)
worst, n_cmp, have_b = 0.0, 0, False
for p in sorted((B / NEW / SUB).glob('firnice_firnledger_*.dat')):
    L = pd.read_csv(p, sep=r'\s+', comment='#', header=None)
    have_b = L.shape[1] >= 7
    if not have_b:
        break
    L.columns = ['elev', 'year', 'we', 'z', 'n', 'old', 'b', 'lim'][:L.shape[1]]
    for e, s in L.groupby('elev'):
        s = s.sort_values('year')
        y0, b = int(s.year.iloc[0]), s.b.to_numpy()
        bref = (s.we.iloc[0] - b[0]) / fl.AGE if s.n.iloc[0] >= fl.AGE else 0.0   # steady state inferred from year 1
        out = fl.run(b, y0, b_ref=bref)
        diff = np.abs(np.array([o[0] for o in out]) - s.we.to_numpy())
        worst = max(worst, float(diff.max())); n_cmp += len(diff)
ok['3 ledger = Python reference'] = have_b and worst < 5e-3
print(f'3: {n_cmp} band-years compared, max |IDL - Python| = {worst:.5f} m w.e.' if have_b else '3: no balance column in the log (older log format): not checked')
# 4: wet depth (year Y) against the deepest limit the model used in Y (ledger of Y-1 + the month's snow)
over, tot = 0, 0
for p in sorted((B / NEW / SUB).glob('firnice_firnledger_*.dat')):
    g = p.stem.split('_')[-1]
    L = pd.read_csv(p, sep=r'\s+', comment='#', header=None, names=['elev', 'year', 'we', 'z', 'n', 'old', 'b', 'lim'])
    W = pd.read_csv(B / NEW / SUB / f'firnice_wetdepth_{g}.dat', sep=r'\s+', comment='#', header=None, names=['elev', 'year', 'wd', 'lost'])
    m = W.merge(L, on=['elev', 'year'])
    over += int((m.wd > m.lim + 1e-6).sum()); tot += len(m)
ok['4 wet depth within the limit used'] = over == 0
print(f'4: {tot} firn band-years; wet depth beyond the logged limit: {over}')
for g, lo, hi in (('02103', 3295, 3365), ('01112', 3275, 3345), ('02596', 3415, 3535)):
    L = pd.read_csv(B / NEW / SUB / f'firnice_firnledger_{g}.dat', sep=r'\s+', comment='#', header=None).iloc[:, :4]
    L.columns = ['elev', 'year', 'we', 'z']
    z = L[L.elev.between(lo, hi)].assign(dec=lambda x: (x.year // 10) * 10).groupby('dec').agg(we=('we', 'mean'), z=('z', 'mean'))
    W = pd.read_csv(B / NEW / SUB / f'firnice_wetdepth_{g}.dat', sep=r'\s+', comment='#', header=None, names=['elev', 'year', 'wd', 'lost'])
    w = W[W.elev.between(lo, hi)].assign(dec=lambda x: (x.year // 10) * 10).groupby('dec').wd.mean()
    print(f'5: {g} {lo}-{hi} m by decade: ledger m w.e. / depth m / wet depth m'); print(pd.concat([z, w], axis=1).round(1).T.to_string())
print('GATE:', {k: ('PASS' if v else 'FAIL') for k, v in ok.items()})
