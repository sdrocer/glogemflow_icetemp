"""Vectorised Python port of firn_change.pro (Beer et al. 2026, TC 20).

Arrays are kept in file orientation (row 0 = top line of the ESRI grid).
IDL linear index -1 (last element, used when WHERE returns -1) is the
top-right cell in IDL's rotated layout = file-orientation [0, nc-1].
"""
import os
import sys
import numpy as np

NOVAL = -9999.0
MBDIR = '/scratch_net/vierzack06_fourth/mhuss/MassBalance_all/'


def read_agr(fn):
    with open(fn) as f:
        head = [f.readline() for _ in range(6)]
        vals = np.array(f.read().split(), dtype=float)
    nc = int(float(head[0][6:46]))
    nr = int(float(head[1][6:46]))
    xll = float(head[2][10:50])
    yll = float(head[3][10:50])
    cs = float(head[4][9:49])
    return vals[:nc * nr].reshape(nr, nc), dict(nc=nc, nr=nr, xll=xll, yll=yll, cs=cs, head=head)


def run_glacier(glacier, tran=(1970, 2025), time_firn=20, mbdir=MBDIR):
    y0, y1 = tran
    L = y1 - y0 + 1
    b_ref, hdr = read_agr(f'{mbdir}{glacier}/massbalance/results/massbalmeas_{y1:04d}.grid')
    nr, nc = b_ref.shape
    glacier_mask = b_ref != NOVAL
    bal = np.full((L, nr, nc), NOVAL)
    present = np.zeros(L, bool)
    for i in range(L):
        fn = f'{mbdir}{glacier}/massbalance/results/massbalmeas_{y0 + i:04d}.grid'
        if os.path.exists(fn):
            bal[i] = read_agr(fn)[0]
            present[i] = True

    firn = np.zeros((nr, nc))
    ageacc = np.full((L, nr, nc), NOVAL)
    firnyear = ageacc.copy()
    firntot = ageacc.copy()
    firn_duration = np.zeros((nr, nc), int)
    year_firn_lost = np.full((nr, nc), NOVAL)
    firn_area_ts = np.zeros(L, int)
    dur_hist = np.zeros((L, nr, nc), int)
    simple_hist = np.zeros((L, nr, nc))

    for i in range(L):
        yr = y0 + i
        b = bal[i]
        acc = b > 0
        abl = (b <= 0) & (b != NOVAL) & (firn > 0)
        # simple cumulative tracker (gates the layer removal)
        firn[acc] += b[acc]
        firn[abl] += b[abl]
        firn[~glacier_mask] = NOVAL
        neg = (firn < 0) & (firn > -100)
        if neg.any():
            firn[neg] = 0
        else:
            firn[0, nc - 1] = 0  # IDL: firn(-1)=0 when WHERE finds nothing
        # new layer for this year
        c = np.full((nr, nc), NOVAL)
        c[acc] = b[acc]
        firnyear[i] = c
        e = np.full((nr, nc), NOVAL)
        e[acc] = yr
        ageacc[i] = e
        # emptying: full |b| compared against each older layer, newest first
        if abl.any():
            active = abl.copy()
            absb = np.abs(b)
            for j in range(i - 1, -1, -1):
                if not active.any():
                    break
                fy = firnyear[j]
                gt = absb > fy
                rem = active & gt
                red = active & ~gt
                fy[rem] = 0
                ageacc[j][rem] = NOVAL
                fy[red] = fy[red] + b[red]
                active = rem
        # layers older than time_firn become ice
        old = (ageacc < yr - time_firn) & (ageacc != NOVAL)
        firnyear[old] = 0
        ageacc[old] = NOVAL
        # age-aware total
        prevtot = firntot[i - 1].copy() if i > 0 else np.zeros((nr, nc))
        fyi = firnyear[:i + 1]
        valid = fyi != NOVAL
        cnt = valid.sum(0)
        s = np.where(valid, fyi, 0.0).sum(0)
        upd = (firn >= 0) & (cnt > 0)
        firntot[i][upd] = s[upd]
        curtot = firntot[i]
        # duration / disappearance
        ll = (prevtot > 0) & (curtot <= 0)
        year_firn_lost[ll & (firn_duration >= 2)] = yr
        firn_duration[ll] = 0
        firn_duration[curtot > 0] += 1
        firn_area_ts[i] = (firn_duration >= 2).sum()
        dur_hist[i] = firn_duration
        simple_hist[i] = firn

    firn_final = firntot[L - 1].copy()
    firn_final[~glacier_mask] = NOVAL
    age = np.full((nr, nc), NOVAL)
    pos = firn_final > 0
    a = ageacc[:, pos]
    v = a != NOVAL
    n = v.sum(0)
    mage = np.where(v, y1 - a, 0).sum(0) / np.maximum(n, 1)
    tmp = age[pos]
    tmp[n > 0] = mage[n > 0]
    age[pos] = tmp
    tsf = np.full((nr, nc), NOVAL)
    ll = (year_firn_lost != NOVAL) & (firn_final <= 0)
    tsf[ll] = y1 - year_firn_lost[ll]
    return dict(firntot=firntot, firn_final=firn_final, age=age, time_since_firn=tsf,
                firn_area_ts=firn_area_ts, dur=dur_hist, simple=simple_hist, bal=bal,
                present=present, glacier_mask=glacier_mask, hdr=hdr,
                years=np.arange(y0, y1 + 1))


if __name__ == '__main__':
    out = sys.argv[1]
    for g in sys.argv[2:]:
        r = run_glacier(g)
        np.savez_compressed(os.path.join(out, f'port_{g}.npz'), **{k: v for k, v in r.items() if k != 'hdr'},
                            hdr=np.array([r['hdr'][k] for k in ('nc', 'nr', 'xll', 'yll', 'cs')]))
        fc = (r['firn_final'] > 0).sum() * 100. / (r['firn_final'] >= 0).sum()
        print(g, 'firn cover %d %%' % int(fc), 'years present', r['present'].sum())
