"""Decision 1 of the firn-ledger plan: the firn model's own layer removal against a mass-conserving one, on the
GLAMOS glaciers. The original (validated port) wipes every layer thinner than |b| without reducing |b| and then
takes the full |b| from the next thicker one; the mass-conserving variant removes exactly |b|, newest first."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = os.environ.get('FIRNMODEL_DIR', '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel')   # port_out/ in, mc/ out
import firn_change_port as fp
NOVAL = fp.NOVAL
src = open(fp.__file__).read()
start = src.index('        # emptying: full |b|'); stop = src.index('        # layers older than time_firn become ice')
mc = '''        # emptying, mass-conserving: exactly |b| removed, newest layers first
        if abl.any():
            rem = np.where(abl, np.abs(b), 0.0)
            for j in range(i - 1, -1, -1):
                if not (rem > 0).any():
                    break
                fy = firnyear[j]
                has = (fy != NOVAL) & (fy > 0) & (rem > 0)
                take = np.where(has, np.minimum(rem, np.where(has, fy, 0.0)), 0.0)
                fy[has] = fy[has] - take[has]
                rem = rem - take
                gone = has & (fy <= 1e-12)
                fy[gone] = 0
                ageacc[j][gone] = NOVAL
'''
ns = {}
exec(compile(src[:start] + mc + src[stop:].split("if __name__ == '__main__':")[0], 'mc_port', 'exec'), ns)
rows = []
for g in sys.argv[1:]:
    q = np.load(f'{S}/port_out/port_{g}.npz')
    m = ns['run_glacier'](g)
    cs = float(q['hdr'][4])
    for lab, r in (('original', q), ('mass_conserving', m)):
        ft, dur = r['firntot'], r['dur']
        for k, y in enumerate(r['years']):
            f = ft[k]; on = f >= 0
            rows.append(dict(glacier=g, rule=lab, year=int(y), firn_km2=float(((dur[k] >= 2) & on).sum() * cs * cs / 1e6),
                             mean_firn_we=float(np.mean(f[on & (f > 0)])) if (on & (f > 0)).any() else 0.0,
                             firn_mass_we_km2=float(np.sum(np.where(on & (f > 0), f, 0.0)) * cs * cs / 1e6)))
    print(g, 'done', flush=True)
D = pd.DataFrame(rows)
D.to_csv(f'{S}/mc/removal_compare_annual.csv', index=False)
D['decade'] = (D.year // 10) * 10
T = D.pivot_table(index=['glacier', 'decade'], columns='rule', values=['firn_km2', 'mean_firn_we'], aggfunc='mean').round(3)
T.to_csv(f'{S}/mc/removal_compare_decades.csv')
print(T.to_string())
