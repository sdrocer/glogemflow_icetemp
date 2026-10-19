"""Check 5 of the firn-ledger plan (plan of 2026-10-19, step 1): GloGEM's firn ledger against the firn model
(firn_change.pro, Beer et al. 2026) on the GLAMOS glaciers, taken apart one step at a time.

  A  the firn model: GLAMOS cells, its own rules (original removal, removal only where its running sum is
     positive, firn counted after 2 years in a row, zero firn in 1970)
  B  the same cells with the ledger's rules, one at a time: B1 mass-conserving removal, B2 no running-sum gate,
     B3 firn counted from the first year; B = all three
  C  the ledger's rules on the GLAMOS band balance (one series per RGI7 10 m band): the effect of averaging the
     balance over the band before the rules; C' the same, started from the steady state of the first ten
     observed years (as GloGEM starts from b_ref)
  D  the ledger GloGEM keeps (firnice_firnledger_<id>.dat of a FRHL run): the same rules, fed GloGEM's balance
Step A-B: the rules; B-C': the aggregation and the start; C'-D: the input (GloGEM's mass balance).

Check that the comparison itself is right: rules() set to the firn model's, run on the GLAMOS cells, must equal
the port (firntot and the firn duration) exactly.

Measures per glacier and year, on the RGI7 bands with at least 3 GLAMOS cells: the firn line (lowest band where
at least half the area has firn) and the firn [m w.e.] in the borehole bands.

Usage: python ledger_vs_firnmodel.py [OUT_DIR]   (FIRNMODEL_DIR: port_out/ and band_tables_v2/, as band_compare.py)
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from firn_change_port import read_agr, NOVAL, MBDIR  # noqa: E402

FM = os.environ.get('FIRNMODEL_DIR', '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel')
RUN = ('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/bayescal/firn_FRHL_r0.85_i1.017_a1.00_w0/'
       'training_runs/monthly/CentralEurope/PAST/firnice_temperature')
GL = {'corvatsch': '02103', 'hohsaas': '01112', 'aletsch': '02596', 'alphubel': '01150', 'tortin': '01039',
      'sexrouge': '01309', 'fiescher': '02611', 'gorner': '01225', 'morteratsch': '02216'}
FULL = ('hohsaas', 'corvatsch', 'aletsch', 'alphubel')     # GLAMOS grids over the whole run period
BH = {'hohsaas': (3275, 3345), 'corvatsch': (3295, 3345), 'aletsch': (3415, 3635), 'alphubel': (3755, 3825)}
MINCELLS = 3
AGE = 20


def rules(bal, y0, removal='orig', gate=True, two_year=True, start=None, quirk=None, colmask=None, age=AGE):
    """Firn bookkeeping on columns (cells or bands). bal (L, N), NOVAL = no glacier / no data that year.
    removal 'orig' (firn_change.pro) or 'mc' (mass-conserving); gate: removal only where the running sum of the
    balance is positive (firn_change.pro); two_year: firn counted after 2 years in a row; start: per-column
    steady-state balance (age + 1 layers before y0) or None (zero); quirk: flat index the port sets to 0 when no
    running sum is negative (IDL firn(-1) = 0); colmask: columns inside the outline (the port keeps firn only inside
    its 2025 outline). Returns firntot (L, N; NOVAL where never firn), present (L, N)."""
    L, N = bal.shape
    S0 = age + 1 if start is not None else 0
    fy = np.full((S0 + L, N), NOVAL)
    ay = np.full((S0 + L, N), NOVAL)
    if start is not None:
        pos = start > 0
        for k in range(S0):
            fy[k, pos] = start[pos]
            ay[k, pos] = y0 - S0 + k          # oldest first
    inside = np.ones(N, bool) if colmask is None else colmask
    run_sum = np.where(inside, 0.0, NOVAL)
    firntot = np.full((L, N), NOVAL)
    dur = np.zeros(N, int)
    present = np.zeros((L, N), bool)
    for i in range(L):
        yr, j0 = y0 + i, S0 + i
        b = bal[i]
        valid = b != NOVAL
        acc = valid & (b > 0)
        abl = valid & (b <= 0) & ((run_sum > 0) if gate else True)
        if gate:
            run_sum[acc] += b[acc]
            run_sum[abl] += b[abl]
            run_sum[~inside] = NOVAL
            neg = (run_sum < 0) & (run_sum > -100)
            if neg.any():
                run_sum[neg] = 0
            elif quirk is not None:
                run_sum[quirk] = 0
        fy[j0, acc] = b[acc]
        ay[j0, acc] = yr
        if abl.any():
            if removal == 'orig':
                active = abl.copy()
                absb = np.abs(b)
                for j in range(j0 - 1, -1, -1):
                    if not active.any():
                        break
                    gt = absb > fy[j]
                    rem = active & gt
                    red = active & ~gt
                    fy[j, rem] = 0
                    ay[j, rem] = NOVAL
                    fy[j, red] = fy[j, red] + b[red]
                    active = rem
            else:
                rem = np.where(abl, np.abs(b), 0.0)
                for j in range(j0 - 1, -1, -1):
                    if not (rem > 0).any():
                        break
                    has = (fy[j] != NOVAL) & (fy[j] > 0) & (rem > 0)
                    take = np.where(has, np.minimum(rem, np.where(has, fy[j], 0.0)), 0.0)
                    fy[j, has] -= take[has]
                    rem = rem - take
                    gone = has & (fy[j] <= 1e-12)
                    fy[j, gone] = 0
                    ay[j, gone] = NOVAL
        old = (ay < yr - age) & (ay != NOVAL)
        fy[old] = 0
        ay[old] = NOVAL
        v = fy[:j0 + 1] != NOVAL
        cnt = v.sum(0)
        s = np.where(v, fy[:j0 + 1], 0.0).sum(0)
        upd = (cnt > 0) & inside & ((run_sum >= 0) if gate else True)
        prev = firntot[i - 1] if i > 0 else np.zeros(N)
        firntot[i][upd] = s[upd]
        cur = firntot[i]
        dur[(prev > 0) & (cur <= 0)] = 0
        dur[cur > 0] += 1
        present[i] = (dur >= 2) if two_year else (cur > 0)
    return firntot, present


def band_ids(g, years, mask, shape):
    """RGI7 band centre of every GLAMOS cell and year, from the nearest DEM (as band_compare.py)."""
    ref, h = read_agr(f'{MBDIR}{g}/massbalance/results/massbalmeas_2025.grid')
    dems = []
    for f in sorted(glob.glob(f'{MBDIR}{g}/topo/dhm_*.grid')):
        z, hz = read_agr(f)
        if z.shape == ref.shape and abs(hz['xll'] - h['xll']) < 1 and abs(hz['yll'] - h['yll']) < 1:
            dems.append((int(re.findall(r'(\d{4})', os.path.basename(f))[-1]), z))
    mids = np.empty((len(years),) + shape)
    for i, y in enumerate(years):
        _, z = min(dems, key=lambda d: (abs(d[0] - y), -d[0]))
        mids[i] = np.floor(z / 10.) * 10 + 5
    return mids


def to_bands(firntot, present, mids, mask, years):
    """Per year and band: mean firn over all cells [m w.e.] and the share of cells with firn."""
    rows = []
    for i, y in enumerate(years):
        m = mids[i][mask]
        f = np.where(firntot[i] > 0, firntot[i], 0.0)[mask]
        p = present[i][mask]
        d = pd.DataFrame(dict(band_mid=m, f=f, p=p)).groupby('band_mid').agg(n=('f', 'size'), we=('f', 'mean'),
                                                                             share=('p', 'mean'))
        d['year'] = y
        rows.append(d.reset_index())
    return pd.concat(rows, ignore_index=True)


def firn_line(d):
    v = d[d.share >= 0.5]
    return v.band_mid.min() if len(v) else np.nan


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(FM, 'check5')
    os.makedirs(out, exist_ok=True)
    lines, thick, exact = [], [], []
    for g, gid in GL.items():
        q = np.load(f'{FM}/port_out/port_{g}.npz')
        years, mask, bal = q['years'], q['glacier_mask'], q['bal']
        nr, nc = mask.shape
        y0 = int(years[0])
        B = bal.reshape(len(years), -1)
        # check: the firn model's rules reproduce the port exactly
        cm = mask.ravel()
        ft, pr = rules(B, y0, quirk=nc - 1, colmask=cm)
        e_ft = np.nanmax(np.abs(np.where(q['firntot'].reshape(len(years), -1) == NOVAL, 0, q['firntot'].reshape(len(years), -1))
                                - np.where(ft == NOVAL, 0, ft)))
        e_pr = int(((q['dur'].reshape(len(years), -1) >= 2) != pr).sum())
        exact.append(dict(glacier=g, max_abs_firn=e_ft, presence_mismatch=e_pr))
        mids = band_ids(g, years, mask, (nr, nc))
        variants = {'A_firn_model': dict(),
                    'B1_mc_removal': dict(removal='mc'),
                    'B2_no_gate': dict(gate=False),
                    'B3_first_year': dict(two_year=False),
                    'B_ledger_rules': dict(removal='mc', gate=False, two_year=False)}
        bands = {}
        for name, kw in variants.items():
            ft, pr = rules(B, y0, quirk=nc - 1 if kw.get('gate', True) else None, colmask=cm, **kw)
            bands[name] = to_bands(ft.reshape(len(years), nr, nc), pr.reshape(len(years), nr, nc), mids, mask, years)
        # C: the ledger's rules on the GLAMOS band balance (band_compare.py's mb_mwe, RGI7 bands)
        t = pd.read_csv(f'{FM}/band_tables_v2/annual_bands_{g}_{gid}.csv')
        t = t[t.glamos_cells >= MINCELLS]
        mb = t.pivot_table(index='year', columns='band_mid', values='mb_mwe').reindex(years)
        bandmids = mb.columns.to_numpy()
        Mb = np.where(np.isfinite(mb.to_numpy()), mb.to_numpy(), NOVAL)
        obs = Mb[:, :]
        first10 = np.array([np.mean([x for x in obs[:, k] if x != NOVAL][:10]) if (obs[:, k] != NOVAL).any() else 0.0
                            for k in range(obs.shape[1])])
        for name, st in (('C_band_zero_start', None), ("C'_band_steady_start", first10)):
            ft, pr = rules(Mb, y0, removal='mc', gate=False, two_year=False, start=st)
            d = pd.DataFrame(dict(year=np.repeat(years, len(bandmids)), band_mid=np.tile(bandmids, len(years)),
                                  we=np.where(ft > 0, ft, 0.0).ravel(), share=pr.ravel().astype(float),
                                  n=np.tile(t.groupby('band_mid').glamos_cells.max().reindex(bandmids).to_numpy(), len(years))))
            bands[name] = d
        # D: GloGEM's ledger (FRHL run), on the same bands
        L = pd.read_csv(f'{RUN}/firnice_firnledger_{gid}.dat', sep=r'\s+', comment='#', header=None,
                        names=['band_mid', 'year', 'we', 'z', 'n_layers', 'oldest', 'b', 'lim'])
        L = L[L.band_mid.isin(bandmids)]
        bands['D_glogem_ledger'] = L.assign(share=(L.we > 0).astype(float))[['year', 'band_mid', 'we', 'share']]
        ok_bands = set(bandmids)
        for name, d in bands.items():
            d = d[d.band_mid.isin(ok_bands)]
            for y, dy in d.groupby('year'):
                lines.append(dict(glacier=g, variant=name, year=int(y), firn_line=firn_line(dy),
                                  present=bool(np.any(mb.loc[y].notna())) if y in mb.index else False))
            if g in BH:
                lo, hi = BH[g]
                z = d[d.band_mid.between(lo, hi)].groupby('year').we.mean()
                for y, v in z.items():
                    thick.append(dict(glacier=g, variant=name, year=int(y), firn_we=v))
        print(g, 'done', flush=True)
    X = pd.DataFrame(exact)
    LN = pd.DataFrame(lines)
    TK = pd.DataFrame(thick)
    X.to_csv(f'{out}/check_exact.csv', index=False)
    LN.to_csv(f'{out}/firnline_by_year.csv', index=False, float_format='%.1f')
    TK.to_csv(f'{out}/borehole_band_firn_by_year.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 250, 'display.max_rows', 400, 'display.max_columns', 30)
    print('check (firn model rules vs the port):'); print(X.to_string(index=False))
    LN = LN[LN.present & LN.year.between(1972, 2025)].copy()
    LN['decade'] = LN.year // 10 * 10
    fl = LN.pivot_table(index=['glacier', 'decade'], columns='variant', values='firn_line', aggfunc='mean')
    fl.to_csv(f'{out}/firnline_by_decade.csv', float_format='%.0f')
    print('\nfirn line by decade [m] (years with GLAMOS data, 1972-2025)'); print(fl.round(0).to_string())
    TK = TK[TK.year.between(1972, 2025)].copy()
    TK['decade'] = TK.year // 10 * 10
    tk = TK.pivot_table(index=['glacier', 'decade'], columns='variant', values='firn_we', aggfunc='mean')
    tk.to_csv(f'{out}/borehole_band_firn_by_decade.csv', float_format='%.2f')
    print('\nfirn in the borehole bands by decade [m w.e.]'); print(tk.round(1).to_string())


if __name__ == '__main__':
    main()
