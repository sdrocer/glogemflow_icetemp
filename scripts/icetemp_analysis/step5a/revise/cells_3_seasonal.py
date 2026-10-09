md(r'''
## 3. The month-to-month error term (T3 only)

**Why T3 needs it.** The model error S treats all rows of one point as fully correlated (they
share a depth). Without another term, a point's months may differ from each other only by the
small measurement error, so every month-to-month wiggle of the misfit counts as strong evidence:
T3's misfit ratio is then 26 instead of about 1 (table below).

**What was measured (Stage 2).** For the 1,503 T3 rows dated to a single month, every pair of
months of the same point: half the squared difference of their misfits (model − borehole), the
*semivariogram* γ. γ = 1 °C² means the two months' misfits typically differ by about 1.4 °C. It
was measured at design 13 and over all 32 runs (they agree within about 10 %), and a covariance K
was fitted to it, pair by pair.

**The term in plain words:**
- **Within one borehole**, the misfit swings with the seasons near the surface. The variance of
  the swing halves every 2.4 m of depth (its size, the standard deviation, every 4.9 m). About a
  quarter of it repeats every year (the same season brings the same error); the rest is forgotten
  within about 7 months.
- **Ice points** get a quarter of firn's variance (ice only conducts). **Firn of entities at or
  above 4,000 m** gets 5.4 times as much: the model melts that cold firn in hot summers, and real
  meltwater soaks in unevenly.
- **Two boreholes at the same point** (all at Monte Rosa) disagree by an extra 0.45 °C², and
  misfits **drift** by about 0.035 °C per year (the model warms faster than the boreholes at
  Colle Gnifetti and Col du Dôme).
- **Nearby depths move together**, over 1.2 m near the surface and 200 m below 20 m.
- **All parts only spread a point's months around the point's own mean.** K adds nothing to the
  uncertainty of a point's mean; S still carries it. One small difference remains (next cell):
  every month row keeps the full measurement error σ² of its borehole entry, independent between
  months, so the mean of a point's T3 rows has a somewhat smaller measurement error than its T2
  row. T1 and T2 themselves are unchanged.
''')

code(r'''
p = ta.SEASONAL
z = np.array([2.0, 5.0, 10.0, 20.0])
display(pd.DataFrame([
    ('variance within one borehole, firn below 4,000 m', f"{p['a']} °C² × exp(−depth ÷ {p['zs']} m)"),
    ('  the same at 2 / 5 / 10 / 20 m depth', ' / '.join(f'{v:.2f}' for v in p['a'] * np.exp(-z / p['zs'])) + ' °C²'),
    ('factor for ice points', f"{p['ice']}"),
    ('factor for firn of entities at band ≥ 4,000 m', f"{p['high']}"),
    ('share of the swing that repeats every year', f"{p['w']}"),
    ('memory of the rest', f"{p['tau']} months"),
    ('extra variance between boreholes of one point', f"{p['u0']} °C²"),
    ('slow drift of the misfit', f"{p['beta']} °C per year (depth-correlation length {p['ld']:.0f} m)"),
    ('depth over which months move together', ' / '.join(f'{l} m' for l in p['lens'])
     + ' at ' + ' / '.join(f'{d:g}' for d in p['lz'][:-1]) + f" / ≥ {p['lz'][-1]:g} m depth"),
], columns=['part of K (ta.SEASONAL)', 'value']).set_index('part of K (ta.SEASONAL)'))
''')

code(r'''
r3_ = TG['T3'].rows
w_ = r3_['count'].to_numpy(float)
w_ = w_ / pd.Series(w_).groupby(r3_.t2.to_numpy()).transform('sum').to_numpy()
s2_ = r3_.sigma.to_numpy(float) ** 2
d3 = pd.Series(w_ ** 2 * s2_).groupby(r3_.t2.to_numpy()).sum()      # measurement variance of a point's mean, T3
d2 = pd.Series(s2_).groupby(r3_.t2.to_numpy()).first()             # the same for its T2 row
print(f"measurement variance of a point's mean, T3 / T2: median {np.median(d3 / d2):.2f}, smallest {np.min(d3 / d2):.3f}; "
      f"with S added (D + S), T3 / T2: smallest {np.min((d3 + ta.S_MODEL_VARIANCE) / (d2 + ta.S_MODEL_VARIANCE)):.3f}; "
      f"largest sigma² {s2_.max():.2f} °C²")
V13 = np.diag(ta.date_spread_var(r3_, obs, s2['Y3'][ta.BEST_RUN]))
r13 = TG['T3'].Y[13] - TG['T3'].y
shr = {}
for share in (0.0, 0.5, 1.0):
    C_ = ta.error_cov(r3_, seasonal=ta.seasonal_cov, shared_d=share) + V13
    bad = []
    for e_ in np.unique(r3_.e):
        j = np.flatnonzero(r3_.e.to_numpy() == e_)
        try:
            np.linalg.cholesky(C_[np.ix_(j, j)])
        except np.linalg.LinAlgError:
            bad.append(e_)
    ok_ = ~np.isin(r3_.e.to_numpy(), bad)
    m_ = ta.misfit_ratio(r13[ok_], C_[np.ix_(ok_, ok_)], r3_[ok_].reset_index(drop=True))
    shr[f'{share:.0%} shared'] = {'entities Sigma cannot be factorised for': ', '.join(r3_.entity[r3_.e == b].iloc[0] for b in bad) or '–',
                                  'ratio at design 13': m_['ratio'], 'point means': m_['level'], 'month-to-month': m_['month_to_month']}
print("T3 with part of each borehole's measurement error shared between a point's months (design 13, D + S + K + V); "
      f"K at 30 m depth is only {ta.SEASONAL['a'] * np.exp(-30 / ta.SEASONAL['zs']):.4f} °C²")
display(pd.DataFrame(shr).T.round(3))
got('3.shared50_ratio', shr['50% shared']['ratio at design 13'])
''')

md(r'''
The measurement error of a point's mean is the same in T3 as in T2 for most points (median ratio
1.0) and down to 0.05 times as large for points read in many months; with the model error S added
it is at least 0.95 times as large, because σ² ≤ 0.81 °C² is small against S = 3.5 °C². With half
of each borehole's measurement error shared between a point's months, T3's misfit ratio goes from
0.66 to 0.73 (its month-to-month part from 0.68 to 0.80). Sharing all of it is not an option: K
almost vanishes at depth, so nothing but the measurement error separates two months of a deep
point; fully shared, the month-to-month part explodes (about 20,000) and one entity's Σ cannot be
factorised at all. So the measurement error stays independent between months, as before; section
5.1 shows what the half-shared version does to T3's cost.
''')

code(r'''
rep2 = json.loads((STEP5A / 'stage2' / 'report.json').read_text())
wb = rep2['within_borehole']
lags = ['1', '2', '3', '4-5', '6', '7-8', '9-11', '12', '13-23', '24-59', '60-119']
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.7), sharey=True)
for ax, band in zip(axes, ['0-5 m', '5-15 m']):
    npts = wb[f"('points', '{band}')"]
    use = [l for l in lags if (npts.get(l) or 0) >= 5]             # lag classes with at least 5 points
    x = np.arange(len(use))
    ax.plot(x, [wb[f"('d13', '{band}')"][l] for l in use], 'o-', color=INK, ms=5, label='measured, design 13')
    ax.plot(x, [wb[f"('fit', '{band}')"][l] for l in use], 's-', color=MUTED, ms=4.5,
            label='the term K, plus measurement error')
    ax.set_xticks(x, use)
    ax.set_xlabel('months between the two readings')
    ax.set_title(f'readings {band} deep')
    ax.set_ylim(0, 2.45)
axes[0].set_ylabel('γ  [°C²]')
axes[0].legend(loc='upper right')
fig.suptitle('Figure 1. How much the misfit of one borehole changes between two months', x=0.01, ha='left', fontsize=10.5)
fig.tight_layout()
plt.show()
''')

md(r'''
Near the surface the misfit follows a yearly cycle: it differs most between months half a year
apart and less between the same months of different years (γ falls from 1.2 °C² at 7–8 months
to 0.5 °C² at 12 months). The term follows these swings only roughly: at 0–5 m it is too high at
3 months (2.13 against 1.09 °C² measured) and at 12 months (1.34 against 0.54 °C²), so it does
not reproduce the dip at 12 months. At 5–15 m it is too low at 3 months (0.46 against 1.96 °C²):
the June and September 1994 Grenz profiles differ there by summer meltwater that soaked in
unevenly. Below 30 m the misfit barely changes between months
(γ ≤ 0.003 °C² within 2 years, less than glenglat's own σ²). The lag classes with fewer than
5 points are not drawn (0–5 m and 5–15 m at 60–119 months); all classes are in Stage 2's
`report.log`.
''')

code(r'''
t1, t2, t3 = TG['T1'], TG['T2'], TG['T3']
rmse2 = np.sqrt(np.mean((t2.Y - t2.y) ** 2, axis=1))          # T2 RMSE of each of the 32 runs
order = np.argsort(rmse2)
best5_rmse, med_run = [int(i) for i in order[:5]], int(order[len(order) // 2])
t3m = t3lit.rows[~t3lit.rows.range_only.to_numpy()].reset_index(drop=True)   # rows dated to one month
Y3m, y3m = t3lit.Y[:, ~t3lit.rows.range_only.to_numpy()], t3lit.y[~t3lit.rows.range_only.to_numpy()]
spike = np.load(STEP5A / 'stage2' / 'spike_rows_month.npy')                  # six Colle Gnifetti melt rows
spike_keys = set(zip(t3m.t2.to_numpy()[spike], t3m.ym.to_numpy()[spike]))
sp3 = np.array([(a, b) in spike_keys for a, b in zip(t3.rows.t2, t3.rows.ym_first)]) & (t3.rows.kind == 'month').to_numpy()

sets = {  # name: rows, model (runs x rows), observed, D+S, K, add the date spread V
    'T2': (t2.rows, t2.Y, t2.y, ta.error_cov(t2.rows), None, False),
    'T3 (T3r rows)': (t3.rows, t3.Y, t3.y, ta.error_cov(t3.rows), ta.seasonal_cov(t3.rows), True),
    'T3 as first specified': (t3lit.rows, t3lit.Y, t3lit.y, ta.error_cov(t3lit.rows), ta.seasonal_cov(t3lit.rows), False),
    'T3 month rows only': (t3m, Y3m, y3m, ta.error_cov(t3m), ta.seasonal_cov(t3m), False)}


def ratio_at(name, run, with_k=True):
    rows, Y, y, DS, K, add_v = sets[name]
    C = DS + (K if (with_k and K is not None) else 0.0)
    if with_k and add_v:
        C = C + np.diag(ta.date_spread_var(rows, obs, s2['Y3'][run]))       # the run's own spread, as Stage 2
    return ta.misfit_ratio(Y[run] - y, C, rows)


tab = []
C1 = t1.Sigma
r1 = lambda run: float((t1.Y[run] - t1.y) @ np.linalg.solve(C1, t1.Y[run] - t1.y) / len(t1.y))
tab.append(dict(target='T1', rows=len(t1.y), **{'D + S only': r1(13), 'with K': r1(13)},
                **{'level part': np.nan, 'month-to-month part': np.nan},
                **{'best five runs': f"{min(r1(i) for i in best5_rmse):.2f}–{max(r1(i) for i in best5_rmse):.2f}",
                   f'median run ({med_run})': r1(med_run)}))
for name in sets:
    m0, m13 = ratio_at(name, 13, with_k=False), ratio_at(name, 13)
    b5 = [ratio_at(name, i)['ratio'] for i in best5_rmse]
    tab.append(dict(target=name, rows=m13['n'], **{'D + S only': m0['ratio'], 'with K': m13['ratio'],
                    'level part': m13['level'], 'month-to-month part': m13['month_to_month'],
                    'best five runs': f'{min(b5):.2f}–{max(b5):.2f}', f'median run ({med_run})': ratio_at(name, med_run)['ratio']}))
ratios = pd.DataFrame(tab).set_index('target')
print(f'misfit ratio r^T Sigma^-1 r / n at design 13 (best five runs by T2 RMSE: {best5_rmse}; median run {med_run})')
display(ratios)
C3 = t3.Sigma                                        # D + S + K + V(run 13)
r3 = t3.Y[13] - t3.y
no_spike = ta.misfit_ratio(r3[~sp3], C3[np.ix_(~sp3, ~sp3)], t3.rows[~sp3].reset_index(drop=True))
print(f"T3 without the six Colle Gnifetti melt rows ({sp3.sum()} rows): ratio {no_spike['ratio']:.3f}")
got('s2.T3r_ratio_noK', ratios.loc['T3 (T3r rows)', 'D + S only']); got('s2.T3r_ratio', ratios.loc['T3 (T3r rows)', 'with K'])
got('s2.T3r_level', ratios.loc['T3 (T3r rows)', 'level part']); got('s2.T3r_m2m', ratios.loc['T3 (T3r rows)', 'month-to-month part'])
got('s2.T3lit_ratio_noK', ratios.loc['T3 as first specified', 'D + S only']); got('s2.T3lit_ratio', ratios.loc['T3 as first specified', 'with K'])
got('s2.T3m_ratio', ratios.loc['T3 month rows only', 'with K']); got('s2.T1_ratio', ratios.loc['T1', 'with K'])
got('s2.T2_ratio', ratios.loc['T2', 'with K']); got('s2.T3r_median_run', ratios.loc['T3 (T3r rows)', f'median run ({med_run})'])
got('s2.T3r_no_spike', no_spike['ratio'])
''')

code(r'''
by = []
zb = ta.depth_band(t3.rows.depth)
for lab, mask in ([(f'depth {b}', zb == b) for b in bands]
                  + [(s, (t3.rows.system_short == s).to_numpy()) for s in sorted(t3.rows.system_short.unique())]):
    m = ta.misfit_ratio(r3[mask], C3[np.ix_(mask, mask)], t3.rows[mask].reset_index(drop=True))
    by.append(dict(group=lab, rows=m['n'], points=m['points'], ratio=m['ratio'], level=m['level'],
                   month_to_month=m['month_to_month']))
    got(f's2.by_{lab}', m['ratio'])
sens = []
for lab, change in [('as proposed', {}), ('no high-firn factor', dict(high=1.0)), ('no between-borehole part', dict(u0=0.0)),
                    ('no drift', dict(beta=0.0)), ('yearly-cycle share 0', dict(w=0.0)), ('yearly-cycle share 0.5', dict(w=0.5)),
                    ('depth lengths × 0.5', dict(lens=tuple(0.5 * np.array(p['lens'])))),
                    ('depth lengths × 2', dict(lens=tuple(2 * np.array(p['lens'])))),
                    ('variance × 0.5', dict(a=p['a'] * 0.5)), ('variance × 2', dict(a=p['a'] * 2))]:
    Cv = ta.error_cov(t3.rows) + ta.seasonal_cov(t3.rows, change) + np.diag(ta.date_spread_var(t3.rows, obs, s2['Y3'][13]))
    m = ta.misfit_ratio(r3, Cv, t3.rows)
    sens.append(dict(variant=lab, ratio=m['ratio'], month_to_month=m['month_to_month']))
    got(f's2.sens_{lab}', m['ratio'])
print('T3 at design 13: each group on its own, then with K changed')
display(pd.DataFrame(by).set_index('group').round(3))
display(pd.DataFrame(sens).set_index('variant').round(3))
display(pd.DataFrame(rep2['spike_rows']).rename(columns={'obs': 'borehole [°C]', 'model_d13': 'model, design 13 [°C]'})
        .set_index('entity'))
''')

md(r'''
**Reading.**
- **With K, T3 passes the honesty check:** its misfit ratio at design 13 is 0.663 (level part
  0.648, month-to-month part 0.676) and stays inside 0.5–2 at the five best runs by T2 RMSE
  (0.63–0.73) and at the median run (0.62). T1 and T2 do not change (1.99 and 0.635).
- **T3 as first specified fails it** (0.356): its repeated range-dated readings turn the model's own
  month-to-month changes into "misfits". That is why T3 means the T3r rows.
- **The term is honest on average but uneven.** It is right near the surface (0–5 m: 1.41), where
  the seasonal fingerprint lives; generous deeper (0.36–0.62) and at Tête Rousse, Hohsaas,
  Corvatsch and Bossons (0.14–0.53), so their month contrasts get less weight than they could;
  tight at Hintereisferner (1.11). Most of this is real site-to-site difference. The term is also
  the same in every season, while the misfits are not (Stage 2's season-pair table, pairs at most
  12 months apart): it over-covers winter pairs (0–5 m: 1.47 against 0.25 °C² measured; 5–15 m:
  0.51 against 0.09) and summer-to-autumn pairs at 0–5 m (1.28 against 0.59), but under-covers
  summer-to-autumn pairs at 5–15 m (0.27 against 1.31 °C²), where the meltwater arrives.
- **Six rows of Colle Gnifetti borehole CG05-1** (table above) are left out of the fit of K but
  kept in every ratio: in the hot summers of 2008 and 2015 the model melts this cold firn; the
  borehole shows no such warming. Without them the ratio is 0.625.
- **The parts are what matters most:** K is not very sensitive to its own details (0.52–0.89 in
  most variants), but dropping the high-firn factor or the between-borehole part doubles the
  month-to-month ratio (1.34 and 1.50).
''')
