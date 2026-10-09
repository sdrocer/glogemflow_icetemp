md(r'''
## 7. The decision rule, row by row

Every number below is computed in this notebook (surface, cubic surface, real runs, fitted error
model) or read from Stage 3 (emulator, quadratic surface); f_adv = 1.0, T3 = the T3r rows. A check
passes only if every method agrees; **not resolved** means the methods fall on both sides of the
limit. The rule itself is unchanged; two rows are added that it did not foresee: the emulator's E
varying with the parameters (section 4.1) and the error model fitted to the misfits (section 5.4,
the plan's sensitivity run).
''')

code(r'''
emu_d = {r['target']: r['d_main'] for r in S3B['consistency']['variants'] if r['method'] == 'EMU' and r['E'] == 'E'}
ho_rs = {n: HO['RS'][n][0] for n in TG}
ho_em = {n: hor[f'EMU/E/{n}'] for n in TG}
ho_q = {n: hor[f'RSQ3/E/{n}'] for n in TG}
ho_run = {n: HO['RUNS'][n][0] for n in TG}
FIT_LABS = [lab for lab in EMT.index if lab.startswith('fitted')]
WORD = {True: 'PASS', False: 'FAIL', None: 'NOT RESOLVED'}


def verdict(values, limit=0.10):
    """True if every estimate is within the limit, False if none is, None if they straddle it."""
    v = np.asarray(list(values), float)
    return True if (v <= limit).all() else False if (v > limit).all() else None


V, txt = {}, {}
for n in TG:
    s, e = ID2.loc[(n, 'surface')], ID2.loc[(n, 'emulator')]
    V[(n, 'told apart')] = bool(s['told apart'] and e['told apart'])
    txt[(n, 'told apart')] = (f"{WORD[V[(n, 'told apart')]]}: correlation {s['correlation']:.2f} / {e['correlation']:.2f}, "
                              f"f_ins spread {s['spread f_ins']:.2f} / {e['spread f_ins']:.2f}; runs {s['runs passing']} / {e['runs passing']}")
    g_ = gtab.loc[(n, 'emulator, E varying')]
    txt[(n, 'told apart, E varying')] = (f"{'pass' if g_['told apart'] else 'fail'}: correlation {g_['correlation']:.2f}, "
                                         f"f_ins spread {g_['spread f_ins']:.2f}")
    if n == 'T1':
        V[(n, 'consistent')] = None
        txt[(n, 'consistent')] = 'does not apply (level only)'
        V[(n, 'RMSE, August error model')] = V[(n, 'RMSE, fitted error model')] = None
        txt[(n, 'RMSE, August error model')] = (f"reference: {ho_rs['T1']:.3f} °C (emulator {ho_em['T1']:.3f}, "
                                               f"real runs {ho_run['T1']:.3f})")
        txt[(n, 'RMSE, fitted error model')] = (f"reference: {EMT.loc[FIT_LABS, 'held-out T1'].min():.3f}–"
                                               f"{EMT.loc[FIT_LABS, 'held-out T1'].max():.3f} °C")
    else:
        d = dist[f'{n}: level vs shape']
        V[(n, 'consistent')] = bool(d <= 3 and emu_d[n] <= 3)
        txt[(n, 'consistent')] = f"{WORD[V[(n, 'consistent')]]}: {d:.2f} sd apart (emulator {emu_d[n]:.2f})"
        costs = {'cubic': ho_rs[n] - ho_rs['T1'], 'emulator': ho_em[n] - ho_em['T1'], 'quadratic': ho_q[n] - ho_q['T1'],
                 'real runs': ho_run[n] - ho_run['T1']}
        V[(n, 'RMSE, August error model')] = verdict(list(costs.values()) + list(carried[n]))
        txt[(n, 'RMSE, August error model')] = (
            f"{WORD[V[(n, 'RMSE, August error model')]]}: cubic {costs['cubic']:+.3f}, emulator {costs['emulator']:+.3f}, "
            f"quadratic {costs['quadratic']:+.3f}, real runs {costs['real runs']:+.3f} °C; with the stand-in's error near "
            f"the fit {carried[n].min():+.2f} to {carried[n].max():+.2f} °C")
        fc = EMT.loc[FIT_LABS, f'{n} − T1']
        V[(n, 'RMSE, fitted error model')] = verdict(fc)
        txt[(n, 'RMSE, fitted error model')] = (f"{WORD[V[(n, 'RMSE, fitted error model')]]}: {fc.min():+.2f} to "
                                               f"{fc.max():+.2f} °C ({len(fc)} fits)")
    h = hon.loc[n]
    V[(n, 'honest')] = bool(0.5 <= h['ratio'] <= 2.0)
    txt[(n, 'honest')] = (f"{WORD[V[(n, 'honest')]]}{', at the limit, by selection' if n == 'T1' else ''}: "
                          f"{h['ratio']:.3f} at run {h['best run']}; fitted optimum {h['ratio there']:.2f}; "
                          f"5 best runs {h['their ratios']} ({h['of them in 0.5–2']} of 5 inside)")
    a_s, a_e = ID3.loc[(n, 'surface')], ID3.loc[(n, 'emulator')]
    V[(n, 'f_adv freed')] = bool(a_s['f_adv passes'] and a_e['f_adv passes'])
    literal = bool(a_s['spread f_adv'] <= 0.5 and a_s['told apart'])
    txt[(n, 'f_adv freed')] = (f"{'YES' if V[(n, 'f_adv freed')] else 'NO'}: f_adv spread {a_s['spread f_adv']:.2f} / "
                               f"{a_e['spread f_adv']:.2f}, runs {a_s['runs with f_adv ≤ 0.5']} / {a_e['runs with f_adv ≤ 0.5']}"
                               + ('; the surface alone at the best point would pass' if literal and not V[(n, 'f_adv freed')] else ''))


def outcome(n, rmse_key):
    if not V[(n, 'told apart')]:
        return '4 (nothing separates f_ref and f_ins)'
    fails = [k for k in ('consistent', 'honest') if V[(n, k)] is False]
    r = V[(n, rmse_key)]
    if r is None:
        return '2 by the stand-ins, 3 by the real runs: not resolved' if fails == ['consistent'] else 'not resolved'
    if r is False:
        return '3 ("RMSE" fails: you decide)'
    return '1 (all pass)' if not fails else '2 (only "consistent" fails)' if fails == ['consistent'] else 'see the rows'


rule = {'told apart': "f_ins spread ≤ 0.5 and correlation ≤ 0.8, surface / emulator (E taken at the surface's best point), "
                      "best point and ≥ 4 of 5 runs",
        'told apart, E varying': "not part of the rule: exact grid fit with the emulator's E varying as in the pipeline (4.1)",
        'consistent': 'level-only and shape-only optima ≤ 3 sd apart',
        'RMSE, August error model': 'held-out RMSE ≤ T1 + 0.10 °C by every method, including the stand-in error near the fit',
        'RMSE, fitted error model': 'the same with an error model fitted to the misfits (5.4; the plan\'s sensitivity run)',
        'honest': 'misfit ratio 0.5–2 at the best run (August error model)',
        'f_adv freed': 'f_adv spread ≤ 0.5 with 3 parameters and f_ref, f_ins still told apart (read like "told apart")'}
DEC = pd.DataFrame({n: {k: txt[(n, k)] for k in rule} for n in TG})
DEC.insert(0, 'passes if', pd.Series(rule))
DEC.loc['outcome, August error model'] = ['', *[outcome(n, 'RMSE, August error model') for n in TG]]
DEC.loc['outcome, fitted error model'] = ['', *[outcome(n, 'RMSE, fitted error model') for n in TG]]
DEC.index.name = 'check'
with pd.option_context('display.max_colwidth', 220):
    display(DEC)
cut = {m: 1 - ID2.loc[('T3', m), 'spread f_ins'] / ID2.loc[('T2', m), 'spread f_ins'] for m in ('surface', 'emulator')}
d32 = {'cubic surface': ho_rs['T3'] - ho_rs['T2'], 'emulator': ho_em['T3'] - ho_em['T2'],
       'quadratic surface': ho_q['T3'] - ho_q['T2'], 'real runs': ho_run['T3'] - ho_run['T2']}
d32_fit = EMT.loc[FIT_LABS, 'T3 − T1'] - EMT.loc[FIT_LABS, 'T2 − T1']
higher = sum(v > 0 for v in d32.values())
print(f"T3 against T2. The rule: T3 only if it cuts the f_ins spread by a further 20 % at no higher RMSE cost.\n"
      f"  f_ins spread cut by {100 * cut['surface']:.0f} % (surface) and {100 * cut['emulator']:.0f} % (emulator)\n"
      f"  held-out RMSE, T3 − T2 (August error model): " + ', '.join(f'{k} {v:+.3f}' for k, v in d32.items())
      + f" °C; higher by {higher} of 4 methods; 90 % interval (cubic) {bs['a_boreholes_refit T3-T2']['p05']:+.2f} to "
      f"{bs['a_boreholes_refit T3-T2']['p95']:+.2f} °C\n"
      f"  with the fitted error model: {d32_fit.min():+.3f} to {d32_fit.max():+.3f} °C\n"
      f"  -> read literally, the rule selects {'T2' if higher else 'T3'}")
got('7.cut_surface', cut['surface'])
got('7.cut_emulator', cut['emulator'])
''')

md(r'''
### Outcome and recommendation

**The outcome depends on the error model.**
- **With the August error model** (the one the plan fixed in advance): **T2 is outcome 2.** It tells
  f_ref and f_ins apart, its held-out RMSE cost is inside your 0.10 °C by every method (+0.07 °C
  by the stand-ins, +0.01 °C by the real runs), and it passes the honesty check; only "consistent"
  fails. **T3 is not resolved between outcome 2 and
  3:** inside the margin by the stand-ins (+0.08 °C), outside it by the real runs (+0.23 °C), and the
  stand-ins' error near its optimum (up to about +0.15 °C) is larger than the margin left.
- **With an error model fitted to the misfits** (section 5.4): **both are outcome 3,** "RMSE
  fails, you decide": +0.20 to +0.37 °C for T2 and +0.19 to +0.29 °C for T3. That error model
  matches the misfits better (its parts pass the honesty check, the August model's do not), and
  lets the shape steer more of the fit, which costs RMSE because shape and level disagree.
- **With T1, f_ins cannot be calibrated** (outcome 4 for T1 alone). **f_adv stays at 1.0** with
  every target.
- **A risk the rule did not foresee:** with the emulator's own uncertainty E varying, as the
  pipeline computes it, neither T2 nor T3 is told apart on this 29-run design (section 4.1).

**Recommendation: build T2, with the August error model and f_adv = 1.0, and document the
level/shape conflict** like Colle Gnifetti's deep ice. At T2's fit the high firn of Monte Rosa
(≥ 3,900 m) is −0.85 °C too cold on average and Col du Dôme −0.58 °C (by the cubic surface), and
the held-out RMSE is expected to be 0.07 °C higher than with T1 (90 % interval −0.12 to +0.25 °C).
Why T2:
- **it is what the rule fixed in advance selects.** T3 cuts the f_ins spread by another 33–34 %
  (the rule asks 20 %), but not "at no higher RMSE cost": T3 − T2 is +0.005 °C by the cubic
  surface, +0.003 °C by the emulator and +0.219 °C by the real runs; only the quadratic surface
  gives T3 the lower cost (−0.058 °C);
- its RMSE check passes by every method; T3's is not resolved;
- at its fit the high firn is less cold than at T3's (Monte Rosa −0.85 against −1.30 °C, Col du
  Dôme −0.58 against −1.01 °C);
- it is the smaller code change (no month rows in the emulator, no month-to-month term).

**T3 stays an option for you, not the rule's choice.** For it: it separates the parameters better
(f_ins spread left 0.12 against 0.19); it stays clearly told apart where T2 sits at the limits
(without the top 3 m T2 reaches correlation 0.77–0.80, without Monte Rosa 0.76) and with E held
fixed (T3 9 of 9 cases, T2 fails 5 of the 6 cases at its optima); with the fitted error model it is
the cheaper of the two (by 0.01–0.09 °C). Against it: its RMSE cost is not resolved and at T3's
best real run (run 23) the high firn is −1.90 °C (Monte Rosa) and −1.66 °C (Col du Dôme); with E
varying it does worse than T2. Only runs near its optimum can settle its cost; the campaign's
second stage provides them, and T2 could be switched to T3 then if you wish.

**Two caveats for either target.**
1. *The error model decides the outcome* (above). The August model passes the honesty check as a
   whole (T2 0.61) but not in its parts: the borehole means are under-covered (2.35) and the shape
   within a borehole over-covered (0.44), so a fit leans on the level.
2. *The emulator's own uncertainty E.* With E varying, neither target is told apart on 29 runs: the
   fit is drawn to the edge of the design, where E is largest. Before adoption, 5b and the campaign
   must handle it: runs near the optimum in the second stage, a cap on E or a tighter coverage
   gate, and a check that the campaign emulator's optimum does not sit on the gate's edge.

### What you decide

Decisions that change the outcome:

| # | decision | recommendation | why |
|---|---|---|---|
| 1 | the error model: the August one (3.5 °C², 15 m), or one fitted to the misfits | the August one, as fixed in the plan; the fitted one stays the sensitivity, and both are re-measured on the campaign's 25 glaciers | outcome 2 for T2 with the August model, outcome 3 for both with the fitted one (+0.19 to +0.37 °C); the August model passes the honesty check as a whole, not in its parts |
| 2 | accept the level/shape conflict and build a depth-resolved target (outcome 2), or keep T1 and calibrate f_ref alone (f_ins = 1.0) | build it | with T1, f_ins cannot be calibrated (spread left 0.66; the data alone 0.92); your "at least two parameters" needs T2 or T3 |
| 3 | T2 or T3 | T2, as the rule selects; T3 only if runs near its optimum confirm its cost | section above |
| 4 | your RMSE margin (0.10 °C) | keep it | T2's expected cost is +0.07 °C, but it exceeds 0.10 °C in 33 % of the resamples; with the fitted error model the margin is exceeded anyway |

Decisions with a clear recommendation (the checks of section 6):

| # | decision | recommendation | why |
|---|---|---|---|
| 5 | readings at the melting point: ordinary, or one-sided | ordinary | the pile-up is gone (15 % of modelled points against 10 % of readings; August: 75 % against 8 %); the one-sided version changes nothing (check c) |
| 6 | Colle Gnifetti below 40 m in or out of the fit | out, reported separately | no effect either way; a known limitation then cannot steer the parameters (check f) |
| 7 | the depth reference under snow | no shift now; check again with the campaign's emulator | the surface moves by at most 0.03 prior sd; the larger moves are the emulator's own (check d) |
| 8 | option C, extra error for readings before 1979 | not needed | old readings are not clearly worse once the Hintereisferner tongue is set aside (check g) |
| 9 | f_adv | fixed at 1.0 | it fails its rule with every target; read literally with the surface alone at the best point, T3's f_adv would pass (section 4) |

The per-glacier offset δ stays out of the fit (your earlier decision).

**For 5b and the campaign:** `explained_variance: 0.999` in the campaign yaml (your file; at 0.99 the
emulator loses T2's f_ins information); the handling of E (caveat 2); runs near the optimum in the
second stage (for "told apart" with E, and for T3's cost if you want T3); the misfit-ratio parts
reported next to the whole.

**Limitations: what may not carry over.** Everything here rests on 7 glaciers, 32 runs (29 for the
stand-ins) and 6 folds. Monte Rosa holds 55 % of the points and most of both the separation and the
conflict. The "told apart" numbers are expectations from slopes and grid fits on stand-ins, not
calibrations. The cost intervals are wide (T2 − T1: −0.12 to +0.25 °C). The campaign will compare
116 entities on 25 glaciers, 18 of which 5a never saw.
''')

md(r'''
## 8. Reproduction check

Each number recomputed above, against the value in the stage reports (the precision of the report
decides what counts as equal: half a unit of its last digit). Five emulator numbers are read from
Stage 3a, not recomputed; for those the check confirms only the transcription. The second table
checks the numbers this revision quotes in its text against the values computed above.
''')

code(r'''
R = {  # key: (reported value, decimals, where)
    's1.rmse_d13_s2': (2.1286, 4, 'Stage 1: design 13 T2 RMSE alps_s2'), 's1.rmse_d13_r7': (2.3185, 4, 'Stage 1: alps_r7'),
    's1.rows_T1': (66, 0, 'Stage 1: T1 rows'), 's1.rows_T2': (771, 0, 'Stage 1: T2 rows'),
    's1.rows_T3lit': (3417, 0, 'Stage 1: T3 rows'), 's1.rows_T3r': (1671, 0, 'Stage 1: T3r rows'),
    's1.range_readings': (168, 0, 'Stage 1: range-dated readings'), 's1.range_rows_lit': (1914, 0, 'Stage 1: rows they fill'),
    's1.rmse_alps_s2_T1': (2.096, 3, 'Stage 1: RMSE d13 s2 T1'), 's1.rmse_alps_s2_T3': (1.913, 3, 'Stage 1: s2 T3'),
    's1.rmse_alps_s2_T3r': (2.404, 3, 'Stage 1: s2 T3r'), 's1.rmse_alps_r7_T1': (2.317, 3, 'Stage 1: r7 T1'),
    's1.rmse_alps_r7_T3': (2.142, 3, 'Stage 1: r7 T3'), 's1.rmse_alps_r7_T3r': (2.743, 3, 'Stage 1: r7 T3r'),
    's2.T1_ratio': (1.99, 2, 'Stage 2: T1 ratio'), 's2.T2_ratio': (0.635, 3, 'Stage 2: T2 ratio'),
    's2.T3r_ratio_noK': (26.32, 2, 'Stage 2: T3r, D + S only'), 's2.T3r_ratio': (0.663, 3, 'Stage 2: T3r with K'),
    's2.T3r_level': (0.648, 3, 'Stage 2: level part'), 's2.T3r_m2m': (0.676, 3, 'Stage 2: month-to-month part'),
    's2.T3lit_ratio_noK': (14.56, 2, 'Stage 2: T3 as specified, D + S'), 's2.T3lit_ratio': (0.356, 3, 'Stage 2: with K'),
    's2.T3m_ratio': (0.686, 3, 'Stage 2: month rows only'), 's2.T3r_median_run': (0.62, 2, 'Stage 2: median run'),
    's2.T3r_no_spike': (0.625, 3, 'Stage 2: without the 6 spike rows'),
    's2.by_depth 0-5 m': (1.41, 2, 'Stage 2: 0-5 m'), 's2.by_depth 5-15 m': (0.62, 2, 'Stage 2: 5-15 m'),
    's2.by_depth 15-30 m': (0.56, 2, 'Stage 2: 15-30 m'), 's2.by_depth > 30 m': (0.36, 2, 'Stage 2: > 30 m'),
    's2.by_Bossons': (0.53, 2, 'Stage 2'), 's2.by_Corvatsch': (0.50, 2, 'Stage 2'), 's2.by_Hintereisferner': (1.11, 2, 'Stage 2'),
    's2.by_Hohsaas': (0.24, 2, 'Stage 2'), 's2.by_Monte Rosa': (0.78, 2, 'Stage 2'), 's2.by_Tête Rousse': (0.14, 2, 'Stage 2'),
    's2.sens_no high-firn factor': (1.02, 2, 'Stage 2'), 's2.sens_no between-borehole part': (1.11, 2, 'Stage 2'),
    's2.sens_no drift': (0.74, 2, 'Stage 2'), 's2.sens_yearly-cycle share 0': (0.69, 2, 'Stage 2'),
    's2.sens_yearly-cycle share 0.5': (0.72, 2, 'Stage 2'), 's2.sens_depth lengths × 0.5': (0.67, 2, 'Stage 2'),
    's2.sens_depth lengths × 2': (0.73, 2, 'Stage 2'), 's2.sens_variance × 0.5': (0.89, 2, 'Stage 2'),
    's2.sens_variance × 2': (0.52, 2, 'Stage 2'),
    '3a.T1_surface_corr': (0.88, 2, 'Stage 3a'), '3a.T1_surface_spread_ins': (0.66, 2, 'Stage 3a'),
    '3a.T2_surface_corr': (0.58, 2, 'Stage 3a'), '3a.T2_surface_spread_ins': (0.188, 3, 'Stage 3a'),
    '3a.T3_surface_corr': (0.34, 2, 'Stage 3a'), '3a.T3_surface_spread_ins': (0.124, 3, 'Stage 3a'),
    '3a.T1_surface_spread_ref': (0.27, 2, 'Stage 3a'), '3a.T2_surface_spread_ref': (0.09, 2, 'Stage 3a'),
    '3a.T3_surface_spread_ref': (0.06, 2, 'Stage 3a'),
    '3a.T1_emulator_corr': (0.87, 2, 'Stage 3a (read)'), '3a.T2_emulator_spread_ins': (0.230, 3, 'Stage 3a (read)'),
    '3a.T3_emulator_spread_ins': (0.155, 3, 'Stage 3a (read)'),
    '3a.T1_surface_spread_adv': (0.91, 2, 'Stage 3a'), '3a.T2_surface_spread_adv': (0.62, 2, 'Stage 3a'),
    '3a.T3_surface_spread_adv': (0.48, 2, 'Stage 3a'), '3a.T2_emulator_spread_adv': (0.75, 2, 'Stage 3a (read)'),
    '3a.T3_emulator_spread_adv': (0.72, 2, 'Stage 3a (read)'),
    '3a.r7_T1_surface_corr': (0.88, 2, 'Stage 3a: alps_r7'), '3a.r7_T1_surface_spread_ins': (0.68, 2, 'Stage 3a: alps_r7'),
    '3a.r7_T2_surface_corr': (0.58, 2, 'Stage 3a: alps_r7'), '3a.r7_T2_surface_spread_ins': (0.19, 2, 'Stage 3a: alps_r7'),
    '3a.r7_T3_surface_corr': (0.33, 2, 'Stage 3a: alps_r7'), '3a.r7_T3_surface_spread_ins': (0.13, 2, 'Stage 3a: alps_r7'),
    '3b.ho_RS_T1': (2.193, 3, 'Stage 3b'), '3b.ho_RS_T2': (2.264, 3, 'Stage 3b'), '3b.ho_RS_T3': (2.269, 3, 'Stage 3b'),
    '3b.ho_RS_T2-T1': (0.072, 3, 'Stage 3b'), '3b.ho_RS_T3-T1': (0.076, 3, 'Stage 3b'), '3b.ho_RS_T3-T2': (0.005, 3, 'Stage 3b'),
    '3b.ho_RUNS_T1': (2.178, 3, 'Stage 3b'), '3b.ho_RUNS_T2': (2.186, 3, 'Stage 3b'), '3b.ho_RUNS_T3': (2.405, 3, 'Stage 3b'),
    '3b.ins_RS_T1_rmse': (2.126, 3, 'Stage 3b'), '3b.ins_RS_T2_rmse': (2.165, 3, 'Stage 3b'), '3b.ins_RS_T3_rmse': (2.250, 3, 'Stage 3b'),
    '3b.dist_T2: level vs shape': (5.41, 2, 'Stage 3b'), '3b.dist_T3: level vs shape': (6.90, 2, 'Stage 3b'),
    '3b.dist_T3: level vs depth shape': (5.31, 2, 'Stage 3b'), '3b.dist_T3: level vs month-to-month': (7.29, 2, 'Stage 3b'),
    '3b.dist_T3: depth shape vs month-to-month': (2.43, 2, 'Stage 3b'),
    '3b.opt_level alone (T1)_f_ref': (0.530, 3, 'Stage 3b'), '3b.opt_level alone (T1)_f_ins': (1.005, 3, 'Stage 3b'),
    '3b.opt_T2 shape alone_f_ref': (0.240, 3, 'Stage 3b'), '3b.opt_T2 shape alone_f_ins': (0.845, 3, 'Stage 3b'),
    '3b.opt_T3 shape alone_f_ref': (0.200, 3, 'Stage 3b'), '3b.opt_T3 shape alone_f_ins': (0.830, 3, 'Stage 3b'),
    '3b.att_T2_MonteRosa_level': (46.6, 1, 'Stage 3b'), '3b.att_T3_MonteRosa_level': (61.5, 1, 'Stage 3b'),
    '3b.att_T2_all_shape': (53.0, 1, 'Stage 3b'), '3b.att_T3_all_shape': (169.3, 1, 'Stage 3b'),
    '3b.honest_T1': (1.992, 3, 'Stage 3b'), '3b.honest_T2': (0.605, 3, 'Stage 3b'), '3b.honest_T3': (0.588, 3, 'Stage 3b'),
    '3b.honest_T2_level': (2.35, 2, 'Stage 3b'), '3b.honest_T3_level': (2.92, 2, 'Stage 3b'),
    '3c.hef_rows_T2': (41, 0, 'Stage 3c'), '3c.hef_rows_T3': (161, 0, 'Stage 3c'),
    '3c.hef_bias_all rows': (-2.49, 2, 'Stage 3c'), '3c.hef_bias_without HEF tongue': (-2.53, 2, 'Stage 3c'),
    '3c.near_T2': (76, 0, 'Stage 3c'), '3c.near_T3': (124, 0, 'Stage 3c'), '3c.colder_T2': (43, 0, 'Stage 3c'),
    '3c.top3_T3_rows_removed': (206, 0, 'Stage 3c'), '3c.cg_rows_T2': (48, 0, 'Stage 3c'), '3c.cg_rows_T3': (100, 0, 'Stage 3c'),
    '3c.old_T2': (65, 0, 'Stage 3c'), '3c.old_T3': (224, 0, 'Stage 3c'),
}
for s, a, b in (('Monte Rosa', 422, 776), ('Bossons', 124, 152), ('Tête Rousse', 106, 144), ('Hintereisferner', 57, 208),
                ('Corvatsch', 41, 301), ('Hohsaas', 21, 90)):
    R[f's1.sys_T2_{s}'], R[f's1.sys_T3r_{s}'] = (a, 0, 'Stage 1: T2 points by system'), (b, 0, 'Stage 1: T3r rows by system')
for b_, v in zip(bands, (93, 211, 180, 287)):
    R[f's1.T2_depth_{b_}'] = (v, 0, 'Stage 1: T2 points by depth')
for lab, v in (('firn', 488), ('ice', 198), ('unknown', 85)):
    R[f's1.T2_label_{lab}'] = (v, 0, 'Stage 1: T2 points by label')
R['s1.T2_mix_months'], R['s1.T2_mix_seasons'] = (403, 0, 'Stage 1: points mixing months'), (341, 0, 'Stage 1: mixing seasons')
folds_rep = {'T1': (1.73, 2.39, 2.55, 0.62, 2.47, 1.20), 'T2': (1.47, 2.41, 2.36, 0.57, 2.65, 1.21), 'T3': (1.46, 2.41, 2.12, 0.56, 2.69, 1.22)}
for n, vals in folds_rep.items():
    for s, v in zip(SYSTEMS, vals):
        R[f'3b.fold_{n}_{s}'] = (v, 2, 'Stage 3b: per fold')


def compare(ref):
    t = pd.DataFrame([dict(number=k, reported=v, recomputed=GOT.get(k, np.nan), source=src,
                           equal=abs(GOT.get(k, np.nan) - v) <= 0.5 * 10 ** -dp + 1e-9) for k, (v, dp, src) in ref.items()])
    t['difference'] = t.recomputed - t.reported
    return t


rep = compare(R)
read = rep.source.str.endswith('(read)')
print(f'{rep.equal.sum()} of {len(rep)} numbers equal the stage reports at their stated precision '
      f'({(~read).sum()} recomputed here, {read.sum()} read from Stage 3a)')
display(rep[~rep.equal].set_index('number').assign(recomputed=lambda d: d.recomputed.map('{:.5f}'.format),
                                                   difference=lambda d: d.difference.map('{:+.5f}'.format)))
print('not computed:', [k for k in R if k not in GOT] or 'none')
Q = {  # numbers quoted in the text of this revision: (quoted value, decimals, where)
    '4.T1_surface_data_only_spread_ins': (0.92, 2, 'section 4 reading, decision 2'),
    '4.1_T1_emuE_corr': (0.87, 2, '4.1'), '4.1_T1_emuE_spread_ins': (0.70, 2, '4.1'),
    '4.1_T2_emuE_corr': (0.95, 2, '4.1, summary'), '4.1_T2_emuE_spread_ins': (0.84, 2, '4.1, summary'),
    '4.1_T3_emuE_corr': (0.90, 2, '4.1, summary'), '4.1_T3_emuE_spread_ins': (1.42, 2, '4.1, summary'),
    '4.1_T2_fixedE_fail': (5, 0, '4.1, section 7'), '4.1_T3_fixedE_fail': (0, 0, '4.1, section 7'),
    '3.shared50_ratio': (0.73, 2, 'section 3'),
    '5.1_carried_T2_min': (0.05, 2, '5.1'), '5.1_carried_T2_max': (0.07, 2, '5.1'),
    '5.1_carried_T3_min': (0.05, 2, '5.1'), '5.1_carried_T3_max': (0.15, 2, '5.1, section 7, summary'),
    '5.1_T3_cost_kfold': (0.069, 3, '5.1'), '5.1_T3_cost_no_drift': (0.175, 3, '5.1'),
    '5.2_level_MR_level alone': (-0.06, 2, '5.2'), '5.2_level_MR_T2 whole': (-0.85, 2, '5.2, section 7'),
    '5.2_level_MR_T3 whole': (-1.30, 2, '5.2, section 7'),
    '5.3_ratio_at_opt_T1': (1.96, 2, '5.3'), '5.3_ratio_at_opt_T2': (0.59, 2, '5.3'), '5.3_ratio_at_opt_T3': (0.58, 2, '5.3'),
    '5.4_cost_T2_fitted at run 13': (0.20, 2, '5.4, section 7, summary'), '5.4_cost_T2_fitted at run 23': (0.36, 2, '5.4'),
    '5.4_cost_T2_fitted at run 13, in each fold': (0.37, 2, '5.4, section 7, summary'),
    '5.4_cost_T3_fitted at run 13': (0.19, 2, '5.4, section 7, summary'), '5.4_cost_T3_fitted at run 23': (0.27, 2, '5.4'),
    '5.4_cost_T3_fitted at run 13, in each fold': (0.29, 2, '5.4, section 7, summary'),
    '6c.at_pmp_d13_T2': (116, 0, 'check c, section 7'), '6c.at_pmp_d13_T3': (248, 0, 'check c'),
    '7.cut_surface': (0.34, 2, 'sections 4 and 7'), '7.cut_emulator': (0.33, 2, 'sections 4 and 7'),
}
qq = compare(Q)
print(f'{qq.equal.sum()} of {len(qq)} numbers quoted in the text of this revision equal the values computed above')
display(qq[~qq.equal].set_index('number'))
print('not computed:', [k for k in Q if k not in GOT] or 'none')
''')

md(r'''
**Differences from the stage reports, and why** (the cell above lists every recomputed number that
does not match; all others match at the stated precision):
- *Hintereisferner left out, T3:* the Stage 3b report's per-fold table gives 2.12 °C; the value is
  2.11495 °C, which rounds to 2.11, as Stage 3b's own log prints. A rounding slip in the report, no
  consequence.
- *Hintereisferner tongue misfit at the T2 fit without the tongue:* the Stage 3c report gives
  −2.53 °C; the value is −2.52483 °C (its log prints −2.525, and that was rounded again). No
  consequence.
- *Read, not recomputed here:* the emulator's numbers (sections 4 and 5), the bootstrap intervals,
  the quadratic-surface variants and Stage 3c's grid fits.
- *Stage 1's matcher check* compares the pipeline's matcher with a re-implementation of its time
  interpolation only; both read the profile at a depth with the same function (`model_at_depth`).
- *Conventions that differ between stages, not errors:* "the top 3 m" is shallower than 3.0 m in
  Stage 3a (51 T2 points) and 3.0 m or shallower in Stage 3c (58 points); the "5 best runs" are the
  best by each target's own misfit in Stage 3a (and here) but by T2 RMSE in Stage 2 and in Stage 3c's
  check (e); Stage 3c uses the campaign prior for f_ref and exact grid fits, Stage 3a the uniform 0–1
  prior and slope-based expectations; for T3 the date spread V comes from each run's own model in
  Stage 2 and from run 13 in Stage 3. The verdicts do not depend on any of these choices.
- *Rounding in Stage 3b's results.json:* its verdict text says T2 = T1 + 0.071 °C because it
  subtracts two rounded RMSEs; the difference itself is 0.0715 °C, reported as +0.072.
''')

md(r'''
## 9. Where the numbers that are read come from

The notebook recomputes everything that takes seconds and reads the rest from the scratch tree
`/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/`. The scripts that wrote those files
are next to them, not in the repository; before this notebook is committed they should be copied
into the repository (with the small result files) or the scratch tree kept. Rebuilding the cache
alone needs only the repository: `PYTHONPATH=src python scripts/icetemp_analysis/target_analysis.py`.

| file (under step5a/) | what it holds | written by (same folder unless stated) |
|---|---|---|
| cache/datahandler.pkl, obs_rows.pkl, model_rows_alps_{s2,r7}.pkl | borehole rows and model rows | `scripts/icetemp_analysis/target_analysis.py` (repository) |
| stage2/report.json | semivariogram tables, spike rows (section 3) | s2z_report.py (with s2lib.py) |
| stage2/spike_rows_month.npy | the six Colle Gnifetti melt rows | s2d_vfit.py |
| identifiability/results.pkl | emulator numbers of section 4, grid fits with E | s3a_main.py (s3lib.py) |
| identifiability/subsets.csv | rows left out (section 4) | s3c_subsets.py |
| heldout/results.json | emulator and quadratic held-out RMSE, bootstrap, consistency variants | h10_results.py, from h03–h07 (hlib.py) |
| heldout/h08_extra.json, h08b_profile.json, h09_cg_deep.json | runs near the fits, level by group, Colle Gnifetti | h08_extra.py, h08b_profile.py, h09_cg_deep.py |
| checks/*.pkl | checks (a)–(g), Stage 3c's reference fits | c01_baseline.py … c07_old_new.py (c3lib.py) |
| revise/kfold.json | K fitted per fold (5.1) | revise/r01_kfold.py |
| revise/emuE_grid.json | grid fits with E varying (4.1) | revise/r04_emuE_grid.py |
| revise/emuE_fixed.json | fits with E held fixed (4.1) | revise/r05_emuE_fixed.py |

`heldout/h03_grid.pkl` (1.26 GB) is not read and can be deleted.
''')
