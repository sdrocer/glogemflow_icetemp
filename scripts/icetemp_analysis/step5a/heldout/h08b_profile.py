"""Where level and shape disagree, in degrees: at the level-only and shape-only optima (response
surface), (a) the entity-level misfit (depth-weighted entity mean of model - obs) by group, and
(b) the within-entity misfit profile (point misfit minus its entity's depth-weighted mean) by group
and depth band."""
import sys, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 100)
obs, s2, common, design = H.load_all()
T = H.build_targets(obs, s2, common)
t1, t2 = T['T1'], T['T2']
keep = design['keep'].to_numpy()
th3 = design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
rs = H.RespSurf('C2Q').fit(th3[keep], t2.Y[keep])
at = {'level only': (0.53, 1.005), 'T2 shape only': (0.24, 0.845), 'T3 shape only': (0.20, 0.83)}
sysn = t2.system
grp = np.where((sysn == 'Monte Rosa') & (t2.rows['band'].to_numpy() >= 3900), 'Monte Rosa firn >= 3900 m',
               np.where(sysn == 'Monte Rosa', 'Monte Rosa < 3900 m', sysn))
grp_ent = pd.Series(grp).groupby(t2.ent).first().to_numpy()
zb = H.ta.depth_band(t2.rows['depth'].to_numpy(float))
A = t1.A2
ent_pos = {e: i for i, e in enumerate(np.unique(t2.ent))}
row_ent = np.array([ent_pos[e] for e in t2.ent])
lev_rows, prof_rows = [], []
for name, th in at.items():
    r = rs.predict(np.array([th]))[0] - t2.y
    lev = A @ r                                  # entity-level misfit
    dev = r - lev[row_ent]                       # within-entity deviation
    for g in sorted(set(grp)):
        m = grp_ent == g
        lev_rows.append(dict(at=name, group=g, entities=int(m.sum()), mean=float(lev[m].mean()), rms=float(np.sqrt((lev[m] ** 2).mean()))))
        for b in ('0-5 m', '5-15 m', '15-30 m', '> 30 m'):
            mm = (grp == g) & (zb == b)
            if mm.sum():
                prof_rows.append(dict(at=name, group=g, band=b, n=int(mm.sum()), dev=float(dev[mm].mean())))
L = pd.DataFrame(lev_rows)
P = pd.DataFrame(prof_rows)
print('== (a) entity-level misfit, model - obs [C]: mean (rms) over the group\'s entities ==')
LL = L.assign(v=L['mean'].round(2).astype(str) + ' (' + L['rms'].round(2).astype(str) + ')').pivot(index='group', columns='at', values='v')
LL['entities'] = L.groupby('group').entities.first()
print(LL.to_string())
print('\n== (b) within-entity misfit profile: point misfit minus its entity mean [C], mean by group and depth band ==')
PP = P.pivot_table(index=['group', 'band'], columns='at', values='dev').round(2)
PP['points'] = P.groupby(['group', 'band']).n.first()
print(PP.to_string())
json.dump(dict(level=L.to_dict('records'), profile=P.to_dict('records')), open(H.OUT / 'h08b_profile.json', 'w'), indent=1)
