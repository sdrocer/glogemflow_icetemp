"""Gate of the firn-ledger batch of 2026-10-19 (plan development_history/plans/2026-10-19_firn_ledger_plan.md).

Runs at f_ref 0.85 (data/bayescal/<tag>_<variant>_r0.85_i1.017_a1.00_w0): gate set 00786, 00774, 03116 (tag g19 on the
rebased firn-ledger branch, g19pre on the branch before 3b-ii, g360after on GloGEMflow with the 360-day year) and
Hohsaas 01112 alone (tag g19h).
  1 defaults: FRH on the ledger branch = FRH on GloGEMflow (thermal output byte-identical)
  2 3b-ii and the review fixes leave 'fixed' unchanged: FRHL g19 = FRHL g19pre
  3 glaciers without observations: FRHLo = FRHL (byte-identical)
  4 the ledger never touches the mass balance: snow logbooks and firn-line logs of FRHL = FRH, of FRHLo = FRHL; for
    FRHLO (observed thermal mask) the snow logbooks and the mass-balance columns of the firn-line log (ELA, AAR)
  5 the IDL ledger = the Python reference (firn_ledger.py), fed the logged balance, from the logged start
  6 observed bands take the GLAMOS balance: the logged balance of every observed band-year equals firnobs_<id>.dat,
    and the start state equals the observed start balance
Reported: FRHLD against FRHL (temperatures at the borehole points); the thermal firn line of FRHLO against FRHL.

Usage: python gate_b19.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import firn_ledger as fl  # noqa: E402

B = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/bayescal')
SUB = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature'
OBS = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/firnobs')
RUN0 = 1972
MB_LOGS = ('firnice_snowmay', 'firnice_snowbh', 'firnice_firnline')


def run(tag, v):
    return B / f'{tag}_{v}_r0.85_i1.017_a1.00_w0'


def files(r):
    return {p.name: p for p in (r / SUB).glob('*.dat')}


def identical(a, b, only=None):
    fa, fb = files(a), files(b)
    names = sorted(set(fa) | set(fb))
    if only:
        names = [n for n in names if n.startswith(only)]
    diff = [n for n in names if n not in fa or n not in fb or fa[n].read_bytes() != fb[n].read_bytes()]
    return len(names), diff


def ledger(r, gid):
    return pd.read_csv(r / SUB / f'firnice_firnledger_{gid}.dat', sep=r'\s+', comment='#', header=None,
                       names=['elev', 'year', 'we', 'z', 'n', 'old', 'b', 'lim'])


def python_check(r, gid):
    """Largest |IDL - Python| over band-years, Python fed the logged balance from the steady state of the start."""
    L, worst, n = ledger(r, gid), 0.0, 0
    starts = {}
    for e, s in L.groupby('elev'):
        s = s.sort_values('year')
        b = s.b.to_numpy()
        bref = (s.we.iloc[0] - b[0]) / fl.AGE if s.n.iloc[0] >= fl.AGE else 0.0   # start state from year 1
        out = fl.run(b, int(s.year.iloc[0]), b_ref=bref)
        d = np.abs(np.array([o[0] for o in out]) - s.we.to_numpy())
        worst, n = max(worst, float(d.max())), n + len(d)
        starts[e] = bref
    return worst, n, starts


def main():
    ok = {}
    n, d = identical(run('g360after', 'FRH'), run('g19', 'FRH'))
    ok['1 defaults identical'] = not d
    print(f'1: FRH ledger branch vs GloGEMflow: {n} files, differ {d[:5]}')
    n, d = identical(run('g19pre', 'FRHL'), run('g19', 'FRHL'))
    ok['2 fixed column unchanged'] = not d
    print(f'2: FRHL before 3b-ii vs now: {n} files, differ {d[:5]}')
    n, d = identical(run('g19', 'FRHL'), run('g19', 'FRHLo'))
    ok['3 no observations, no change'] = not d
    print(f'3: FRHLo vs FRHL without observation files: {n} files, differ {d[:5]}')
    res4 = []
    for a, b in ((run('g19', 'FRH'), run('g19', 'FRHL')), (run('g19h', 'FRHL'), run('g19h', 'FRHLo'))):
        n, d = identical(a, b, only=MB_LOGS)
        res4.append((a.name, b.name, n, d))
    a, b = run('g19h', 'FRHL'), run('g19h', 'FRHLO')
    n, d = identical(a, b, only=MB_LOGS[:2])
    fl_cols = [pd.read_csv(r / SUB / 'firnice_firnline_01112.dat', sep=r'\s+', comment='#', header=None).iloc[:, :3]
               for r in (a, b)]
    if not fl_cols[0].equals(fl_cols[1]):
        d = d + ['firnice_firnline_01112.dat (ELA/AAR)']
    res4.append((a.name, b.name, n + 1, d))
    ok['4 mass balance untouched'] = all(not r[3] for r in res4)
    for r in res4:
        print(f'4: {r[1]} vs {r[0]}: {r[2]} mass-balance logs, differ {r[3][:4]}')
    res5 = {}
    for r, gids in ((run('g19', 'FRHL'), ('00786', '00774', '03116')), (run('g19h', 'FRHL'), ('01112',)),
                    (run('g19h', 'FRHLo'), ('01112',)), (run('g19h', 'FRHLO'), ('01112',))):
        for g in gids:
            res5[(r.name, g)] = python_check(r, g)
    worst = max(v[0] for v in res5.values())
    ok['5 ledger = Python'] = worst < 5e-3
    print(f'5: {sum(v[1] for v in res5.values())} band-years, max |IDL - Python| {worst:.5f} m w.e.')
    O = pd.read_csv(OBS / 'firnobs_01112.dat', sep=r'\s+', comment='#', header=None, names=['year', 'elev', 'b'])
    O0, Oy = O[O.year == 0].set_index('elev').b, O[O.year >= RUN0]
    msgs = []
    for v in ('FRHLo', 'FRHLO'):
        L = ledger(run('g19h', v), '01112')
        m = L.merge(Oy, on=['elev', 'year'], suffixes=('', '_obs'))
        bad = int((np.abs(m.b - m.b_obs) > 5e-5).sum())
        starts = res5[(run('g19h', v).name, '01112')][2]
        s = pd.Series(starts).reindex(O0.index)
        bad0 = int((np.abs(s - O0.clip(lower=0)) > 5e-4)[s.notna()].sum())
        msgs.append(f'{v}: {len(m)} observed band-years, balance not the observed one: {bad}; '
                    f'start not the observed one: {bad0} of {int(s.notna().sum())} bands')
        ok[f'6 {v} takes the observations'] = bad == 0 and bad0 == 0 and len(m) > 0
    print('6: ' + ' | '.join(msgs))
    # reported
    pa = pickle.load(open(run('g19', 'FRHL') / 'training_runs/design_0000.output.pkl', 'rb'))
    pb = pickle.load(open(run('g19', 'FRHLD') / 'training_runs/design_0000.output.pkl', 'rb'))
    dd = []
    for k in pa:
        if pa[k] is not None and pb[k] is not None:
            a, b = (np.asarray(x['T'] if isinstance(x, dict) and 'T' in x else x, float).ravel() for x in (pa[k], pb[k]))
            if a.shape == b.shape:
                dd.append(b - a)
    D = np.concatenate(dd) if dd else np.array([np.nan])
    print(f'reported: FRHLD - FRHL at {D.size} borehole points: mean {np.nanmean(D):+.3f}, max |d| {np.nanmax(np.abs(D)):.3f} C')
    for v in ('FRHL', 'FRHLO'):
        f = run('g19h', v) / SUB / 'firnice_firnline_01112.dat'
        t = pd.read_csv(f, sep=r'\s+', comment='#', header=None, names=['year', 'ela', 'aar', 'firnline', 'a', 'b'])
        t['dec'] = t.year // 10 * 10
        print(f'reported: {v} thermal firn line, Hohsaas, by decade:', t.groupby('dec').firnline.mean().round(0).to_dict())
    print('GATE:', {k: ('PASS' if v else 'FAIL') for k, v in ok.items()})


if __name__ == '__main__':
    main()
