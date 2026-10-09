md(r'''
## 6. The checks (a)–(g)

Seven questions from the plan, each answered on `alps_s2` with f_adv = 1.0, for T2 and T3. They
were run in Stage 3c as **exact fits on a grid** (not the slope-based expectations of section 4),
with the campaign prior for f_ref (uniform 0.1–1, prior sd 0.260), the quadratic surface and Stage
3c's emulator without E (both held at or below the melting point). So their reference fits are
slightly wider than section 4's: T2 correlation 0.63 and f_ins spread 0.22, T3 0.39 and 0.14
(surface). The tables come from `step5a/checks/`; row counts and a few headline numbers are
recomputed here. Where the surface and the emulator disagree, the emulator's own refitting noise
(retrained on the same runs with 0.01 °C added: f_ins moves by 0.11 prior sd for T2 and 0.21 for
T3, `c08_emu_refit_noise.log`) is the yardstick.

*Columns:* f_ref, f_ins = the fit's mean values; **shift** = how far they move against the fit with
all rows, in prior sd; **distance** = the same move in units of the uncertainty of the fit with all
rows (above about 1 it matters).
''')

code(r'''
CK = STEP5A / 'checks'
ck = lambda f: pickle.load(open(CK / f, 'rb'))
REN = {'n': 'rows', 'mean_f_ref': 'f_ref', 'mean_f_ins': 'f_ins', 'corr': 'correlation', 'spread_f_ins': 'spread f_ins',
       'd_f_ref_prior_sd': 'shift f_ref', 'd_f_ins_prior_sd': 'shift f_ins', 'mahal_base': 'distance'}


def check_table(tab, variants=None, extra=(), methods=('surface', 'emulator')):
    t = tab[tab.method.isin(methods)]
    if variants is not None:
        t = t[t.variant.isin(variants)]
    t = t.assign(target=t.target.replace({'T3r': 'T3'})).set_index(['variant', 'target', 'method'])
    return t[[c for c in REN if c in t] + list(extra)].rename(columns=REN).round(3)


rows2, rows3 = TG['T2'].rows, TG['T3'].rows
grp2, grp3 = ta.group_of(rows2), ta.group_of(rows3)
quad = {n: QuadSurface(theta[keep], TG[n].Y[keep]) for n in ('T2', 'T3')}
pmp = {n: pressure_melting_point(TG[n].rows.depth.to_numpy(float)) for n in ('T2', 'T3')}
model_at = lambda n, f_ref, f_ins: np.minimum(quad[n].coef.T @ quad_basis(((np.array([f_ref, f_ins, 1.0]) - quad[n].c) / SD)[None, :])[0], pmp[n])
''')

md(r'''
### (a) Do firn points and ice points want the same f_ins?

f_ins has two jobs: the firn's conductivity, and since Step 2 the winter snow on ice bands. Fitted
separately (f_ref free in both), firn rows and ice rows are compared; then each borehole group's
preferred f_ins with f_ref held at the joint optimum.
''')

code(r'''
a = ck('c02_firn_ice.pkl')
display(check_table(a['table'], ['all rows', 'label firn', 'label ice', 'model ice band, without Hintereisferner tongue']))
z = a['z'][(a['z'].a.isin(['label firn', 'model firn band']))
           & a['z'].b.isin(['label ice', 'model ice band, without Hintereisferner tongue'])].copy()
z['target'] = z.target.replace({'T3r': 'T3'})
z['firn'] = z.f_ins_a.map('{:.2f}'.format) + ' ± ' + z.sd_a.map('{:.2f}'.format)
z['ice'] = z.f_ins_b.map('{:.2f}'.format) + ' ± ' + z.sd_b.map('{:.2f}'.format)
print('f_ins wanted by firn rows (a) and ice rows (b); difference in combined sd')
display(z.set_index(['target', 'a', 'b', 'method'])[['firn', 'ice', 'z']].round(2))
g = a['groups'][(a['groups'].target == 'T2') & (a['groups'].method == 'surface')]
print('T2, surface: each group\'s preferred f_ins with f_ref at the joint optimum (0.333)')
display(g.set_index('group')[['n', 'label', 'f_ins_data_opt', 'at_bound', 'dchi2_over_range', 'bias_at_joint']]
        .rename(columns={'n': 'points', 'f_ins_data_opt': 'preferred f_ins', 'at_bound': 'at a bound',
                         'dchi2_over_range': 'misfit change over the range', 'bias_at_joint': 'misfit at the joint optimum [°C]'}).round(2))
''')

md(r'''
**Answer.** On average one f_ins serves both: firn rows want 0.79 ± 0.08 and ice rows 0.82 ± 0.25
(T2, surface; 0.13 combined sd apart; in T3 1.0 sd). But only because sites that pull opposite ways
cancel. The Hintereisferner tongue (−2.5 °C at the joint optimum) wants the lowest allowed f_ins,
the Grenz tongue (+1.2 °C) the highest; without the Hintereisferner tongue the ice bands want
1.23 ± 0.29 against 0.79 for the firn (1.5–2.0 combined sd). The firn splits the same way: Colle
Gnifetti wants 1.63, the Monte Rosa high firn 0.525. Ice rows alone learn almost nothing in T2 (f_ins
spread 0.74–0.77, f_ref keeps its whole prior spread).
''')

md(r'''
### (b) With and without the Hintereisferner tongue

The tongue (41 points, 161 T3 rows) is the most biased group and no parameter reaches it.
''')

code(r'''
b = ck('c03_hef.pkl')
display(check_table(b['table'], extra=('rmse T2 all', 'rmse HEF tongue')))
hef2 = grp2 == 'Hintereisferner tongue, ice (< 3000 m)'
for lab, say in (('all rows', 'on all rows'), ('without HEF tongue', 'without the tongue')):
    s = b['fits'][(lab, 'T2')]['surface']
    bias = np.mean(model_at('T2', s['map_f_ref'], s['map_f_ins'])[hef2] - TG['T2'].y[hef2])
    print(f"T2 fit {say:19s} (f_ref {s['map_f_ref']:.3f}, f_ins {s['map_f_ins']:.3f}): tongue misfit {bias:+.3f} °C")
    got(f'3c.hef_bias_{lab}', bias)
print(f"design 13 (real run): tongue misfit {np.mean(TG['T2'].Y[13, hef2] - TG['T2'].y[hef2]):+.2f} °C; rows T2 {hef2.sum()}, "
      f"T3 {(grp3 == 'Hintereisferner tongue, ice (< 3000 m)').sum()}")
got('3c.hef_rows_T2', hef2.sum()); got('3c.hef_rows_T3', (grp3 == 'Hintereisferner tongue, ice (< 3000 m)').sum())
''')

md(r'''
**Answer.** The tongue barely steers the fit: without it f_ins moves by +0.10 (surface) and +0.24
(emulator) prior sd, a distance of 0.44–0.72 in T2 and 0.68–0.97 in T3. Telling the parameters
apart is unchanged (T2 correlation 0.63 → 0.65), and so is the RMSE over all points (2.163 →
2.162 °C). The tongue stays about 2.5 °C too cold at every fitted parameter set: no parameter
reaches its misfit.
''')

md(r'''
### (c) Readings at the melting point

**Your August question.** In August, at strong warming, the model sat at the melting point at 75 %
of the points against 8 % of the readings; a point pinned there tells the fit nothing. The advice
then was to re-measure that share once the likelihood compares whole profiles, and to add a
one-sided (censored) treatment only if it stayed large. So the first answer is the share itself.

The second is a one-sided version, as an approximation: a reading at the melting point is a lower
bound (the model cannot be warmer than 0 °C, adjusted for pressure: depth × −0.000668 °C/m,
GloGEM's clamp), so its misfit is set to zero when the model is at or above min(reading, melting
point) and counted as usual when the model is colder. This is a hinge on single misfits, not a full
censored likelihood, and by construction it can forgive at most about 0.1 °C per reading.
''')

code(r'''
near = {n: TG[n].y >= pmp[n] - 0.1 for n in ('T2', 'T3')}
for n in ('T2', 'T3'):
    above = TG[n].y > pmp[n]
    print(f"{n}: {near[n].sum()} of {len(near[n])} rows within 0.1 °C of the melting point or above it "
          f"({above.sum()} above it, by at most {(TG[n].y - pmp[n])[above].max():.3f} °C)")
    got(f'3c.near_{n}', near[n].sum())
display(check_table(ck('c05_melting.pkl')['table']))
s = ck('c05_melting.pkl')['fits'][('ordinary', 'T2')]['surface']
G = model_at('T2', s['map_f_ref'], s['map_f_ins'])
yb = np.minimum(TG['T2'].y, pmp['T2'])
colder = (G < yb - 1e-9) & near['T2']
print(f"at the T2 optimum (surface): the model is colder than the reading at {colder.sum()} of {near['T2'].sum()} such points, "
      f"by {np.mean((G - yb)[colder]):+.2f} °C on average; those are penalised either way")
got('3c.colder_T2', colder.sum())
for n in ('T2', 'T3'):
    at_d13 = TG[n].Y[13] >= pmp[n] - 0.01
    at_all = (TG[n].Y[ik] >= pmp[n] - 0.01).all(axis=0)
    line = (f"{n}: model at the melting point (within 0.01 °C) at {at_d13.sum()} of {len(at_d13)} rows at design 13 "
            f"({100 * at_d13.mean():.0f} %), at {at_all.sum()} in all 29 runs; readings within 0.1 °C of it: "
            f"{near[n].sum()} ({100 * near[n].mean():.0f} %)")
    if n == 'T2':
        line += f"; at the T2 optimum (surface): {(G >= pmp['T2'] - 0.01).sum()} points"
    print(line)
    got(f'6c.at_pmp_d13_{n}', at_d13.sum())
''')

md(r'''
**Answer.**
- **The share, re-measured: the pile-up at the melting point is gone.** At design 13 the model is
  at the melting point at 116 of the 771 T2 points (15 %), at 101 in every one of the 29 runs, and
  at 115 at the T2 optimum (surface); 76 readings (10 %) are within 0.1 °C of it (August: 75 %
  against 8 %).
  T3: 248 of 1,671 rows (15 %) against 124 readings (7 %).
- **The one-sided version changes nothing:** the optimum moves by at most 0.001 prior sd (0.003 in
  units of the fit's uncertainty), and the misfit at the optimum changes by +0.03 (T2) and +0.97
  (T3); it forgives at most 0.086 °C on 24 points (`c09_strengths.log`). The misfit even rises a
  little, because zeroing single misfits inside a correlated error model is not a proper censored
  likelihood. Since the hinge can forgive at most about 0.1 °C per reading, this test could not
  have found a large effect; the share above is what answers your question.
- **The real misfit is elsewhere:** at 43 of the 76 points the model is about 1.0 °C colder than
  the reading, and those count either way.
- Leaving these points out moves f_ins by +0.11 (surface) and +0.34 (emulator) prior sd; with the
  emulator the correlation then reaches 0.802, just past the rule's limit of 0.8 (distance 1.03),
  about three times the emulator's refitting noise.
- **Recommendation (your decision, section 7): keep them as ordinary readings;** a censored term is
  not needed for T2 or T3.
''')

md(r'''
### (d) The depth reference under snow

The model counts depth from the snow surface, the boreholes from the ice surface. Where a shallow
reading was taken while the model's ice band carried snow (the `firnice_snowbh` logbook), reading
the model at depth + snow depth would compare like with like.
''')

code(r'''
d = ck('c04c_snow_fit.pkl')
print(f"rows at 3 m or shallower read (partly) under modelled snow on an ice band: T2 {d['p2_snow'].sum()}, "
      f"T3 {d['r3_snow_shallow'].sum()}")
display(check_table(d['table'], extra=('rmse T2 shallow under snow',)))
''')

md(r'''
**Answer.** 13 of the 58 T2 points at 3 m or shallower, and 61 of the 206 such T3 rows (mostly
Hintereisferner 2805 m, March–September 1972–74), were read under a modelled snow cover of 2.1 m
on the median (0.1–3.3 m; `c04c_snow_fit.log`, line 5). Shifting them by the snow depth changes the
model there by −0.20 °C on average (0.97 °C RMS; `c04b_snow_shift.log`) and makes the tongue's
misfit there slightly worse (−2.60 → −2.80 °C; `c04c_snow_fit.log`, line 6). Dropping or shifting
them moves the surface optimum by at most 0.03 prior sd. The emulator moves more, and differently:
for T2 its f_ins moves by up to 0.24 prior sd, about twice its refitting noise (0.11); for T3 by
0.22–0.24 prior sd, about its noise (0.21), but 1.4 times its fit's uncertainty, and with all rows
shifted its correlation falls from 0.41 to 0.03. These are effects of the emulator that the surface
does not show. **Recommendation (your decision, section 7): no depth shift in the fit now; check it
again with the campaign's emulator.**
''')

md(r'''
### (e) T3 without the top 3 m

Logger depths drift as the snow buries them or melt lowers the surface. Are the parameters still
told apart without the rows at 3 m depth or shallower?
''')

code(r'''
e = ck('c06_top3m_cg.pkl')
display(check_table(e['table'], ['all rows', 'without top 3 m']))
lin = {}
for n in ('T2', 'T3'):
    t, r = TG[n], ID['alps_s2'][(n, '2p')]
    m = t.rows.depth.to_numpy(float) > 3.0
    pts = [('best point', r['best'])] + [(f'run {i}', np.r_[theta[i, :2], 1.0]) for i in r['runs']]
    res = [posterior(r['jac'](q)[m], t.Sigma[np.ix_(m, m)]) for _, q in pts]
    lin[n] = {'rows': int(m.sum()), 'best point: correlation': res[0][0], 'best point: spread f_ins': res[0][1][1],
              'runs: correlation': f"{min(x[0] for x in res[1:]):.2f}–{max(x[0] for x in res[1:]):.2f}",
              'runs: spread f_ins': f"{min(x[1][1] for x in res[1:]):.2f}–{max(x[1][1] for x in res[1:]):.2f}",
              'runs passing': f"{sum(told_apart(x[0], x[1][1]) for x in res[1:])}/5"}
    got(f'3c.top3_{n}_rows_removed', (~m).sum())
print('Expected from the slopes, as in section 4 (surface), without the rows at 3 m or shallower')
display(pd.DataFrame(lin).T)
''')

md(r'''
**Answer.** T3 still tells them apart: correlation 0.40 → 0.49 (surface) and 0.41 → 0.46
(emulator), f_ins spread 0.14 → 0.17 and 0.24 → 0.19, and the optimum moves by −0.04 and −0.22
prior sd. **T2 only just does:** correlation 0.77 and 0.80, f_ins spread 0.49 and 0.48, right at
the rule's limits (0.8 and 0.5); expected from the slopes (table above) it passes at the best point (0.70 / 0.24)
and at its 5 best runs. The reason: f_ref and f_ins act almost alike below 15 m (similarity of their
effects −0.95) and differently near the surface (−0.62 in T2, −0.38 in T3; `c09_strengths.log`).
Stage 3a's version of this check (top 3 m = shallower than 3.0 m, 7 T2 points fewer) gave T2 0.68 /
0.23 and T3 0.41 / 0.15.
''')

md(r'''
### (f) Colle Gnifetti deeper than 40 m

48 points (100 T3 rows) of 4 entities at 40.8–101 m that carry the climate before 1942, which the
runs do not have (an accepted limitation).
''')

code(r'''
cg2 = np.asarray(dict(ta.group_masks(rows2))['  of it deeper than 40 m (accepted)'], bool)
cg3 = np.asarray(dict(ta.group_masks(rows3))['  of it deeper than 40 m (accepted)'], bool)
print(f"rows: T2 {cg2.sum()}, T3 {cg3.sum()}; design 13 too warm there by {np.mean(TG['T2'].Y[13, cg2] - TG['T2'].y[cg2]):+.2f} °C (T2)")
got('3c.cg_rows_T2', cg2.sum()); got('3c.cg_rows_T3', cg3.sum())
display(check_table(e['table'], ['all rows', 'without CG > 40 m'], extra=('rmse T2 all',)))
h9 = json.loads((STEP5A / 'heldout' / 'h09_cg_deep.json').read_text())
for lab in ('all rows', 'without CG > 40 m'):
    h = h9[lab]['heldout']
    print(f"cubic surface, fits {lab:18s}: held-out RMSE T1 / T2 / T3 " + ' / '.join(f"{h[n]['all_points']:.3f}" for n in TG)
          + f" °C; level vs shape {h9[lab]['d_T2']:.2f} (T2), {h9[lab]['d_T3']:.2f} sd (T3)  [Stage 3b]")
''')

md(r'''
**Answer.** No effect. Without them f_ins moves by −0.04 to −0.07 prior sd, the distance stays
≤ 0.37, correlation and spread are unchanged, and the RMSE over all points goes from 2.163 to
2.161 °C. Held out, the cost over T1 becomes +0.067 (T2) and +0.069 °C (T3), and the level/shape
distance 5.39 and 6.90 sd. In these rows f_ref and f_ins act as one (similarity −0.999), so they say
nothing about telling the two apart (`c09_strengths.log`). **Keeping them out of the fit and
reporting them separately costs nothing** (the plan's recommendation).
''')

md(r'''
### (g) Old against new readings

"Old" = readings before 1979: Hintereisferner (1972–76), Col du Dôme 1973–76 (Bossons 4275 and
4285 m) and, in T3 only, 8 rows of Colle Gnifetti 1977 (Grenzgletscher@4455; in T2 those readings
are pooled with later ones, so their 8 points count as "mixed"). If the old entities fit clearly
worse, your 28 September option C (extra model error for pre-1979 entities) comes back.
''')

code(r'''
old2, new2 = (rows2.ym_last < 197901).to_numpy(), (rows2.ym_first >= 197901).to_numpy()
old3 = (rows3.ym_last < 197901).to_numpy()
print(f'T2 points: old {old2.sum()}, new {new2.sum()}, mixed {(~old2 & ~new2).sum()}; T3 rows: old {old3.sum()}, new {(~old3).sum()}')
got('3c.old_T2', old2.sum()); got('3c.old_T3', old3.sum())
gg = ck('c07_old_new.pkl')
res = gg['residuals']
res = res[res['at'].isin(['design 13', 'surface optimum'])].assign(target=lambda x: x.target.replace({'T3r': 'T3'}))
display(res.pivot_table(index=['target', 'set'], columns='at', values=['bias', 'rmse', 'ratio'], sort=False).round(2))
display(check_table(gg['table']))
''')

md(r'''
**Answer.** Old readings are not clearly worse once the Hintereisferner tongue is set aside: at the
T2 optimum the old points have a misfit ratio of 0.55, the new ones 0.54; without the tongue the old
points' bias is −0.62 °C. The sites measured in both periods do not tell one story: Colle Gnifetti
fits better in 1977 than later (RMSE 1.66 against 3.40 °C at the T3 optimum), and at Col du Dôme the
1973 borehole is 4.0 °C warmer than the model while the 1974–76 one is 1.25 °C colder; those two old
entities, 10 m apart in elevation, disagree with each other (`c07_old_new.log`). In T3 the old rows
have a ratio of 1.06 (new 0.47), almost all from the month-to-month part at Hintereisferner 2805 m
(`c11_old_t3r_split.log`). Leaving the old rows out moves f_ins by +0.10 (surface) and +0.27
(emulator) prior sd; option C would need only 0.05 °C² per old entity for T2, which changes
nothing. **Recommendation (your decision, section 7): option C is not needed.**

**For later stages.** Stage 3c's reference fits are below: the surface, its emulator without E, and
its emulator with E computed once at its first optimum and kept ("emulator_E"). The last one is narrow for
T3 only because its optimum lies on the coverage gate's edge (section 4.1); it is not used for any
verdict.
''')

code(r'''
base = ck('baseline.pkl')
display(pd.DataFrame({(n.replace('T3r', 'T3'), m): {k: v for k, v in base[n][m].items()
                                                    if k in ('map_f_ref', 'map_f_ins', 'corr', 'spread_f_ref', 'spread_f_ins', 'E_top_eig')}
                      for n in base for m in base[n]}).T.rename(columns={'map_f_ref': 'f_ref (optimum)', 'map_f_ins': 'f_ins (optimum)',
                                                                        'corr': 'correlation', 'E_top_eig': 'E: largest pattern [°C²]'}).round(3))
''')
