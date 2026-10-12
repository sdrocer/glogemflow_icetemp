# Snow insulation, a fair comparison, and a calibration target that tells two parameters apart

*Revised 2026-09-30 from three investigations, each re-checked by an independent review, and a
final consistency check. The full version with every number, source and disagreement is
`2026-09-30_snow_comparison_target_plan_full.md` (archived next to this file). The firn-extent plan
of 28 Sep is archived in `development_history/plans/2026-09-28_firn_extent_plan.md`.*

## Context — your three comments

**1. "This is just after the spin-up, not after we calibrate?"** Yes. Every number so far comes from
forward runs with one fixed parameter set each (the 32-run ensemble `alps_r7` spreads them over their
ranges; "best run" = design 13, f_ref 0.52, f_ins 1.09, f_adv 1.05). **No calibration has been run.**

```
 one forward run (~11 min, 7 glaciers)
 spin-up to equilibrium (1942–71 climate) → main run Oct 1971 – Sep 2025 → compare with each reading at its date
 calibration (not run yet): an emulator trained on such runs lets a sampler search BETWEEN them
```

A calibration can only move what its parameters reach. On ice sites all three move the temperature
by 0.1–0.3 °C per prior standard deviation, so no calibration can close Hintereisferner's 4–5 °C.
Such misfits are fixed if they are bugs and documented if they are missing physics.

**2. "Check whether the recent changes improved the model before adding physics; focus on the
calibration."** Agreed, and the plan is built that way:
- **no new physics** — only bug fixes (the code, the comparison or the data do not do what they
  should) and measurements;
- **firn memory is not added** (Step 3: with the correct mass balance it would change nothing the
  boreholes or the calibration see) — your call, since you had asked for it;
- **Step 4 measures** what the fixes bring on the same 32 parameter sets and the same readings;
- **Step 5, the calibration target, is the main work.**

What the recent changes brought (each measured when it landed, with different yardsticks, so they
cannot be summed):

| change | kind | measured effect | now |
|---|---|---|---|
| conduction grid, implicit solver, bed node, column to the borehole (22–24 Sep) | bug fixes | 238 of 776 readings had been invisible | in |
| thermal spin-up: 30-year climate, main-run physics, equilibrium (24–28 Sep) | bug fixes + choices you approved | RMSE 3.75 → 2.51 °C¹ | in |
| each reading at its own month; interpolation between nodes (28 Sep) | comparison fixes | RMSE 3.22 → 2.50 °C¹; 702 → 783 readings | in |
| each glacier's own RGI7 mass balance (29 Sep) | input fix | ensemble best 2.37 → 2.34 °C, median 2.92 → 2.66 °C; ELA follows WGMS | in |
| strain heating (unit fix); surface at 0 °C in melt months (29 Sep) | processes | fit worse | off |

¹ measured with the wrong (RGI6-numbered) mass balance: valid as before/after only.

**How much physics does a flowline model need?** Enough to get each band's thermal regime right
(where the firn is, the surface temperature and its seasons, the winter snow on ice, ice flow, the
geothermal flux). Beyond that, each process brings a parameter that ~800 readings on 7 glaciers cannot
pin down. A new process is worth considering only if a misfit (1) is larger than the parameters can
reach, (2) appears on several glaciers, (3) has a cause with its own signature, and (4) pays in a
paired ensemble test like Step 4.

**3. "Fit two parameters to reduce the RMSE; they must be told apart for the chosen target."**
Exactly — and the target decides it (Step 5). Today's target keeps one mean misfit per borehole
entity; with it f_ins cannot be calibrated. Comparing per depth bin (T2), or per depth bin and month
(T3), tells the two apart (expected correlation after a fit 0.57 or 0.34 instead of 0.87). The catch:
the lowest RMSE and the clearest separation want slightly different parameters; on glaciers left out
of the fit that costs 0.05–0.15 °C.

*Words used:* **entity** = one glacier's boreholes in one 10 m band (e.g. Grenzgletscher@2505);
**point** = an entity's readings pooled into a 1 m depth bin (769 on the 7 glaciers); **prior sd** =
the spread allowed before a fit (f_ref 0.29, f_ins 0.33, f_adv 0.24); **held-out** = scored on a
glacier system the fit did not see; **emulator** = a fast statistical stand-in for GloGEM.

## What we found

1. **Correct mass balance, firn line right.** With RGI7 parameters (icetemp `65de865`) the ELAs match
   the calibration and Hintereisferner follows WGMS (offset −68 m, r 0.83, slope 0.97). `alps_r7`:
   best RMSE 2.36 °C on 769 points (2.34 on the 689 shared with the old ensemble), best runs at
   f_ref 0.38–0.54, inside the range.
2. **Two misfits remain at the best run.** Hintereisferner tongue ice **−4.4 to −5.0 °C** (41 points;
   the earlier "−3.5 °C" also averaged two firn sites — chapter 6 will be corrected). Colle Gnifetti
   below 40 m +1.3 to +1.7 °C (climate before 1942; accepted and documented).
3. **Snow on ice barely insulates, because of a solver bug** (Step 2): thin winter snow passes heat
   5–7× too easily. Fixing it removes about half of the Hintereisferner misfit.
4. **The comparison has two bugs, the held-out test leaks, one borehole sits at the wrong
   elevation** (Step 0).
5. **The reason for firn memory is gone** (Step 3).
6. **Whether f_ref and f_ins can be told apart depends on the target** (Step 5).
7. **Deferred (your decision):** Herron–Langway stage 1 is multiplied by the accumulation rate
   (`firnice_perm_depth.pro`), which the 1980 paper does not do.

## The plan

```
 Step 0  compare model and boreholes fairly (Python, no IDL) ──┐ in parallel
 Step 1  snow logbook (IDL, no physics change)               ──┘
 Step 2  snow conduction bug fix + exact tests           one 11-min run + ~1.5 h spin-up check
 Step 3  firn memory: a decision only (recommended: not built)
 Step 4  measure: the 32-run ensemble again, paired with alps_r7    ~70 min on 6 IDL sessions
         ── stop: you see the table ──
 Step 5  calibration target: 5a analysis → rule → 5b code (Python)
 Then    two-stage training campaign (2 parameters) → KO calibration → adopt or not
```

Steps 0–4 use the 7 test glaciers (RGI7 00773 Bossons, 00777 Tête Rousse, 01112 Hohsaas, 01225
Gorner/Grenz, 01230 Monte Rosa east face, 02103 Corvatsch, 03116 Hintereisferner); single runs use
design 13.

### Step 0 — compare model and boreholes fairly (Python only)

*Why first:* every later number is measured with this ruler, and it costs no IDL time.

**0a. Read the model at the reading's date, not at the end of its month.** *Analogy:* a thermometer
read on 10 July is compared with what the model printed on 31 July — three weeks of summer warming
apart. GloGEM writes each monthly profile after the month's last time step
(`firnice_temperature_model.pro:614`, after the loop ending at 575), and today's matcher compares every
reading of the month with it.

```
     June                      July                      August
 ─────┼─────────────────────────┼─────────────────────────┼───▶
      ● model 30 June    reading 10 July ▲       ● model 31 July  ← today's comparison
      └──────────── interpolate: T = T_before + f·(T_after − T_before) ───┘
```

| symbol | meaning |
|---|---|
| T_before, T_after | model at the reading's depth at the end of the previous month / of the reading's month |
| f | fraction of the month elapsed on the reading's day (0 on the 1st, 1 on the last) |

The correction is 1.0 °C RMS in the top 5 m (up to 4.9 °C), 0.3 °C at 5–10 m, ≤ 0.03 °C below 15 m.
Left in, it would pull a calibrated f_ins down by about 0.3 — several times the calibration's own
uncertainty. On alps_r7: best RMSE 2.36 → 2.34 °C; the Hintereisferner tongue looks 0.1–0.3 °C colder
(about −4.9 °C), because its summer readings are then compared with mid-month values.

**0b. Profiles dated only to a year or range** are compared with January today (16 profiles, 168 of
1722 month rows, all deeper than 10.8 m; error ≤ 0.95 °C). Fix: average the model over the range and
add its spread to the reading's error.

**0c. "Leave one glacier out" must leave the whole glacier system out.** The folds were RGI outlines,
and the Colle Gnifetti saddle straddles two of them; the pipeline's own validation even separates
Gornergletscher@2525 from Grenzgletscher@2525, which read the same model column.

```
   RGI7 01225 (Gorner/Grenz)                  │  RGI7 01230 (Monte Rosa east face)
   Grenz@4445   Grenz@4475   Grenz@4485       │  Grenz@4455   Grenz@4465
   └──── Colle Gnifetti saddle: boreholes 2–33 m apart in elevation ─────┘
```

Fix: one fold per glacier system (01225 + 01230 = Monte Rosa), in the notebook and in
`validation.Validator`. The held-out cost of the shape targets drops from +0.33 °C to +0.05–0.15 °C.

**0d. Borehole l9b (ryser2013, id 292) at the wrong elevation — your decision.** glenglat lists
2015 m, estimated from contours; its published position is 65 m from l9a (2545 m, drilled the day
before). 2015 m is below the model's lowest band, so its 10 points are compared with the wrong ice.
Recommendation: set 2545 m as one documented override in the icetemp data handling (and tell the
glenglat maintainers); alternative: exclude it. No IDL run needed.

**Then re-score alps_r7** from its archived profiles — the baseline for everything below; the set of
points is fixed here.

*Safe because* only the comparison changes, not GloGEM. *Could go wrong:* the model is not linear
within a month (it steps every 10 days), but a line between two month ends is still far closer than
the month end; readings in the run's first month keep today's value. *Tests:* a synthetic sine wave is
recovered within one day of phase; a reading on a month's last day gives today's value; no model column
appears in both training and test of any fold. *Commits:* one icetemp commit per item.

### Step 1 — a logbook of the snow the thermal model sees (no physics change)

*Analogy:* before calling a sleeping bag broken, measure its thickness and put one thermometer outside
and one inside.

```
 ice band in January (Hintereisferner tongue, 2615 m)
   node 0 ● surface: T_surf = min(0 °C, air)             logged: T_surf
          │ snow                                          logged: snow (m w.e.), depth from its mass
   node 1 ●
          │ snow, then ice
   node 2 ● first node under the snow                     logged: T_ice
   R = thickness ÷ conductivity from the surface to the first node under the snow
       R_used = as the solver uses it   R_series = what layered snow and ice should give
```

Three small files per glacier: `firnice_snowclim_<id>.dat` (band × month means), `firnice_snowmay_<id>.dat`
(1 May snow each year, for WGMS) and `firnice_snowbh_<id>.dat` (monthly snow at the borehole bands, for
Step 5). The values are collected inside `firnice_temperature_model.pro` next to its existing output,
which is already silent during the spin-up, so `glogem.pro` is not touched.

**Expected** (from a Python copy of today's column): in January 0.30 m w.e. of snow, 1.17 m deep,
under a −14.7 °C surface; the ice at 2 m only ~5 °C warmer than the surface (~10 °C after Step 2);
R_used ≈ 1.3 against R_series ≈ 7.3 m² K W⁻¹. On 1 May the copy carries ~0.7 m w.e.; WGMS has
0.60 m w.e. at 2600–2650 m (five years, 2010–2015).

**What it decides** — not whether Step 2 happens, but what to expect from it: 0.2–0.4 m w.e. on the
tongue in winter → +2 to +3 °C expected; about 0.1 → smaller (the copy is re-run with the logged snow);
no winter snow → no Hintereisferner expectation; more than 2× off WGMS → a mass-balance question,
written down but not changed here (the mass balance is shared with GMIP4).

*The run* (`snowlog_w0`, ~11 min on vierzack03) reuses alps_r7's configuration and borehole list and
sets the parameters with the ensemble's own writer (0.5245 / 1.0893 / 1.054977). *Safe because* the
check proves it: every `temp_ID*.dat` byte-identical to the archived design 13 (65 profiles), firn line
logs identical. If not, we stop and find out why. *Commits:* one GloGEM, one icetemp.

### Step 2 — let heat cross the snow as the model's own equations say (a bug fix)

*Analogy:* 1 m of insulation in front of 1 m of concrete. Heat must pass both; their resistances add,
and the insulation decides. The solver averaged the two conductivities instead, as if the insulation
were half concrete.

**What is wrong.** Temperatures sit on nodes (0, 1, … 9 m, then 14 m …); heat flows through the gap
between two nodes. Three errors combine:
1. **The snow is stored half a layer too high:** the snow filling 0–1 m goes onto node 0, whose
   temperature is set by the air (its heat capacity is lost there too).
2. **Every gap gets the average conductivity of its two nodes** — fine inside smooth firn, wrong where
   snow meets ice. A partly filled gap averages the densities first, then converts.
3. **The snow depth uses 300 kg/m³** while the snow layers get 250 / 300 / 360 kg/m³ (the column holds
   −17 % to +21 % of the model's snow).

```
          what fills the gap         the solver uses now              after the fix
  node 0 ●
         │  1 m snow, 250 kg/m³      (0.15 + 2.13)/2 = 1.14 W/m/K     0.15 W/m/K
  node 1 ●
         │  ice                      2.13                             2.13
  node 2 ●
```

| snow (m w.e.) | 0.1 | 0.2 | 0.3 | 0.45 | 0.6 | 1.0 |
|---|---|---|---|---|---|---|
| heat crosses the snow too easily by | 4.8× | 7.2× | 5.8× | 3.7× | 1.8× | 1.4× |

**The fix.** In a gap holding snow, resistances add (the textbook rule for layers, Patankar 1980):

  Δz / k_gap = s / k_snow + (Δz − s) / k_col

| symbol | meaning |
|---|---|
| Δz | thickness of the gap between two nodes (m) |
| s | snow inside that gap (m) |
| k_snow | Calonne (2011) conductivity × f_ins, never above ice (W m⁻¹ K⁻¹) |
| k_col | the column below the snow: Yen (1981) ice, ≈ 2.1 |
| k_gap | what the solver uses for that gap |

Plus three consistency repairs: snow depth from its mass and its layers' densities; the heat capacity
of what each node's slice actually holds; snow nodes rising with the ice, not 3.6× faster. **Stays the
same:** firn bands (snow stays inert there, `d072e0f`), no meltwater heat in ice (`97dfbf2`), the
drift spin-up (it sums the solver's own coefficients), the mass balance, T_off.

**Why it is a bug, not physics.** A column in steady state has an exact answer. For 0.15–0.45 m w.e.
of snow today's scheme misses it by 1.1–2.2 °C; the fixed scheme by 2·10⁻⁹ °C. The snow conductivity
was introduced (`1aef7f6`) precisely so that snow would insulate; the solver never carried it out.

**What the column copy shows** (Hintereisferner 2615 m, 1942–71 climate, the model's own snow):

| | now | after the fix | fine-grid reference | no snow |
|---|---|---|---|---|
| annual mean, 10 m | −6.3 °C | −3.3 °C | −3.2 °C | −5.3 °C |
| January, ice at 1 m | −11.0 °C | −5.8 °C | −5.8 °C | |
| effect of f_ins at 10 m (per 0.35) | 0.04 °C | 0.52 °C | | |

Today the snow even *cools* the tongue by 1 °C. The warming holds under other assumptions (ice
emergence +2.7 °C; half / 1.5× the snow +1.2 / +2.7 °C).

**About half of the Hintereisferner misfit.** At the readings' own months and depths: +2.2 to +3.0 °C
against −4.7 °C (about −4.9 after Step 0) — expected afterwards about −2.4 °C. The rest most likely is
melt-season heat (wet snow at 0 °C would add ~1.8 °C) — new physics, not in this plan.

**Who else changes:** every borehole on a band the model had as ice, with snow, before its reading —
251 of the 769 points (Hintereisferner tongue 41, Grenz tongue 52, Gorner tongue 46, Tête Rousse 91,
Hohsaas 21). The other 518 must not move. Rough pooled RMSE: 2.36 → 2.16 °C if only Hintereisferner
moves, up to ~2.3 °C if the Grenz tongue and Tête Rousse warm too (Tête Rousse may get worse).

**Why it matters for Step 5:** f_ins now reaches the ice sites, where f_ref does nothing — a
fingerprint of its own. But the ice stays ~2 °C too cold (asks for a *lower* f_ins) while the firn is
slightly too warm (asks for a *higher* one): Step 5 checks whether both want the same f_ins.

**Tests.** (a) *The IDL code itself:* a one-band IDL test column against the exact steady answer
(within 0.0001 °C; today 1–2 °C off, which must show first) and against the Python copy over a real
year (0.001 °C). (b) *7 glaciers* (`snowfix_w0` vs `snowlog_w0`): the table above site by site; firn
bands within 0.01 °C; stop checks — if Hintereisferner warms by < +1 or > ~+3.5 °C, or the Grenz
tongue by > +2 °C, we compare with the logbook and test column before going on. (c) *Spin-up:* drift
correction vs 3000 plain cycles on 00777, 01225, 01230 + 03116 (~1.5 h), within 0.01 °C for 01225/01230;
if it fails, Step 4 waits. (Tête Rousse alone would test nothing: all its bands are firn in the spin-up.)

*Could go wrong:* the fix can overshoot (Grenz tongue, Tête Rousse too warm). It still stays — a bug
fix — after ruling out an edit error; the result goes to Step 5. Three simplifications start to
matter and are recorded, not changed: start-of-month snow held all month; leftover snow removed only
in January; depth measured from the snow surface in the model but from the ice surface in boreholes.
*Commits:* one GloGEM (the fix), one icetemp (test column, Python copy).

### Step 3 — firn memory: checked; recommendation: leave it out for now

- **The reason is gone.** With the correct mass balance the Grenz tongue is never in the firn area and
  warms, as its boreholes do (20–110 m, 1979–80 → 2008–10: +0.36 to +0.55 °C; with the old parameters
  it cooled by 1–2.6 °C).
- **A store would matter only where a band becomes firn for the first time:** at the 1977, 1978 and
  1981 updates (106 bands, +3.1 °C at 14 m within two years), as the spin-up's firn area moves down to
  the 1970s ELAs. No borehole sits there. The 234 later entries are bands that had been firn before;
  a store would leave them as today.
- **The calibration would not notice:** Tête Rousse and Hohsaas, the only glaciers with boreholes on
  bands moving in and out of the firn area, carry 0.8 % of the information on f_ref and 0.2 % on f_ins.
- **Emptying the store when the mass balance says ice would be wrong:** in 95 cases a band firn for ≥ 10
  years dropped out for one or two years only — firn that old cannot melt that fast.
- **Tête Rousse's warm bias has another cause** (100 % firn from the spin-up until ~2000, meltwater to
  32–67 m): even the best case for a store moves 3155 m only from +1.95 to +1.2 °C.

**Recommendation:** no store now; record this in chapter 6 and revisit it with the Herron–Langway
question, which sets how big a store would be. The snowline calibration (3c of 28 Sep) goes with it:
its premise came from the wrong mass balance. Both reverse earlier wishes of yours, so they need your
OK. (If you want the store anyway: a minimal version is in the full plan — a per-band "savings account"
of firn, kept open while the band is on the glacier; ~15 lines, one run; it will not move the fit.)

### Step 4 — measure what the changes bring, before anything else

*Analogy:* weigh the same 32 samples on the old scale and the repaired one; only the model changed.

```
 design point     1    2    3   …   13   …   32     (the alps_r7 design, seed 11)
 alps_r7          ●    ●    ●        ●        ●  ─┐ both scored with the Step-0 comparison,
 alps_s2 (Step 2) ●    ●    ●        ●        ●  ─┘ on the same points → difference = the fix
```

`alps_s2_w0…w5`: the 32 runs with Step 2 in (~70 min, 6 IDL sessions on vierzack03, profiles archived);
the GloGEM tree unchanged while they run. Reported in one table: pooled RMSE and bias (best run and
median), change per design point, per borehole group, by firn/ice and depth, the strength of f_ref /
f_ins / f_adv on firn and ice points, and where the best runs sit.

**Keep only what pays** (fixed before we look):

| kind of change | stays if |
|---|---|
| bug fix | always; if the fit gets worse, something else was compensating — we look for it, we do not revert |
| artefact fix (e.g. the firn store) | no borehole group worse by > 0.05 °C RMSE |
| new process or parameter | pooled RMSE better by ≥ 0.05 °C (best and median), no group worse by > 0.1 °C, f_ref and f_ins not harder to tell apart |

*Expected (estimates):* Hintereisferner ice about −4.9 → −2.4 °C; pooled RMSE 2.34 → roughly
2.15–2.3 °C (possibly worse if Tête Rousse and the Grenz tongue warm a lot); f_ins on the ice
0.1 → 0.5 °C per prior sd. **Then a stop: you see the table before Step 5's code is written.**

### Step 5 — the calibration target (the main work)

*All Step-5 numbers are estimates from the 32 alps_r7 runs (7 glaciers, 769 points) through a smooth
fit around the best run; 5a repeats everything on the Step-4 ensemble before any code.*

**The idea.** Think of the firn as a room heated in summer and cooled in winter. **f_ref is the
heater** (meltwater heat refreezing in the firn), **f_ins is the wall** (how easily heat moves through
firn and, after Step 2, through the snow on ice). A warmer room can come from either, so the *mean*
temperature cannot tell them apart. Their fingerprints in the profile differ:

```
 depth    heater (f_ref)                            wall (f_ins)
  0–5 m   summer heat added from the top            how far the seasonal swing reaches, and how late
  5–20 m  a warm bulge down to where water stops    less conductive firn damps it shallower, delays it more
  > 20 m  only the long-term mean                   only the long-term mean
 (on ice bands only the wall acts, through the winter snow)
```

**Three targets** (ways to compare model and boreholes):

```
 one entity   depth  Jun-98 Nov-98 Jun-99
                2 m    x      x      x      T1 entity mean     1 number per entity         level
                5 m    x      x      x      T2 depth           1 number per 1 m bin        + shape
               10 m    x      x      x      T3 depth + month   1 number per bin and month  + season
```

| (f_adv fixed at 1.0) | T1 (today) | T2 depth | T3 depth + month |
|---|---|---|---|
| numbers compared | 66 | 769 | 1674 |
| f_ref ↔ f_ins correlation after a fit | 0.87 | 0.57 | 0.34 |
| spread left, f_ref / f_ins (1 = nothing learned) | 0.26 / 0.72 | 0.11 / 0.28 | 0.09 / 0.22 |
| … without Monte Rosa (01225 + 01230) | 0.38 / 0.89 | 0.25 / 0.57 | 0.17 / 0.36 |

- **Today's target cannot calibrate f_ins** (the 0.72 is the prior speaking). Depth-resolved targets
  can, also with the emulator's own uncertainty included and without the top 3 m.
- The separation comes from the **top 15 m of the firn sites** (Monte Rosa above 3900 m, Col du Dôme).
  Month-to-month changes carry a fingerprint of their own.
- **f_adv** keeps 62–94 % of its spread in every target; it stays at 1.0 unless 5a shows otherwise.

**The catch: level and shape disagree.** Fitting the level alone and the shape alone gives optima
6.5 standard deviations apart; the disagreement sits on Monte Rosa's high firn (41 % of the points).
Its price is small: held-out RMSE 2.34 °C with T1, 2.45 with T2, 2.40 with T3. An error model that
absorbs more misfit, or today's per-glacier offset δ inside the fit, makes it costlier (2.56–3.12 °C).
In plain words: for the lowest RMSE the level must steer, to tell the parameters apart the shape must
steer. Steps 2–4 act on exactly this near-surface shape, so we decide on the Step-4 ensemble.

**5a — analysis** (notebook `06_calibration_target.ipynb`, no package code): the table above, the
level/shape split, held-out RMSE per glacier system with a 90 % interval, where the information comes
from; plus: do firn and ice points want the same f_ins; with and without the Hintereisferner tongue;
readings at the melting point (74 points, one-sided or not — your August question); the depth reference
under snow; Colle Gnifetti below 40 m in or out; the three-parameter version (f_adv); old against new
readings (your option C of 28 Sep).

**Decision rule, fixed before we look:**

| check | passes if |
|---|---|
| told apart | f_ins spread left ≤ 0.5 and correlation ≤ 0.8, by two methods, at the best point and 4 of the 5 best runs |
| consistent | level-only and shape-only fits within 3 standard deviations |
| RMSE | held-out RMSE ≤ T1's + 0.10 °C (the margin is yours) |
| honest error model | misfit about as large as the error model expects (ratio 0.5–2) |

Outcomes: **all pass** → build T2 (smaller change; T3 only if it cuts the f_ins spread by another 20 % at
no higher cost). **Only "consistent" fails** → the parameters partly compensate a model error on
Monte Rosa's firn; build if within your margin, document the conflict (your call). **"RMSE" fails** →
you decide: two parameters with the measured cost, or f_ref alone. **Nothing separates them** → f_ref
alone (against your "at least two"), and test a parameter with its own fingerprint (a geothermal-flux
scale, or the percolation depth after the Herron–Langway fix) in a small ensemble first. *On alps_r7
today it would end in outcome 2 with T3, if T3's seasonal error term passes "honest" once re-measured.*

**5b — code** (only for the chosen target; one switch `target = entity | depth | depth_month`,
`entity` reproducing today bit for bit). The likelihood compares every row instead of one mean per
entity:

  log L(θ) = −½ [ rᵀ Σ⁻¹ r + log det Σ ],   r = y − G(θ),   Σ = D + S + E(θ)

| symbol | meaning |
|---|---|
| θ | the two parameters (f_ref, f_ins); f_adv = 1.0 |
| y, G(θ) | measured and emulated temperature of each row (a depth bin, or a bin and month) |
| D | measurement error (glenglat's, ×3 for "estimated" profiles, plus date spread) |
| S | model error, correlated with depth in a borehole (and between months for T3) |
| E(θ) | emulator uncertainty (keeping 99.9 % of the variance) |
| Σ⁻¹, log det Σ | weigh each misfit by the error expected there; count correlated rows about once |

Also: the fixed-parameter bookkeeping (today correct only when the fixed parameter is the last one);
a design over the free parameters only (f_ref 0.1–1.0); validation through GloGEM's emulator instead of
the analytical surrogate, one glacier system per fold; `explained_variance: 0.999` in the campaign
yaml (yours to edit — at 99 % half the f_ins information is lost). *Could go wrong:* 7 glaciers and 32
runs (the ranking T1 < T2 < T3 never changes, absolute numbers do); Monte Rosa dominates (everything is
reported per glacier and without it); Tête Rousse sits at the melting point unless f_ref ≈ 0 (runs
below 0.1 are left out).

### Then — the training campaign and the calibration

```
 f_ins ▲  stage 1: 40 runs over the whole box        stage 2: 20 runs around the stage-1 result
   2.0 │ ·    ·   ·    ·    ·                      │
   1.0 │ ·   ·    ·   ·    ·   ·                   │          · ◉ ·
   0.3 │ ·    ·    ·   ·    ·                      │
       └────────────────────▶ f_ref 0.1–1.0        └────────────────────▶ f_ref
```

- f_ref 0.1–1.0, f_ins 0.3–2.0, f_adv = 1.0; all 116 CentralEurope entities (25 glaciers); spin-up on;
  profiles archived; ~3.5 h + ~2 h on 6 sessions (estimate; one run is timed first, and the runner's
  30-min limit raised).
- Two stages because a single ~40-run design fails the existing "≥ 5 runs close to the result" check in
  73–99 % of draws.
- **Adopt only if** the chains agree (R-hat ≤ 1.05) and match a grid search, the misfit ratio is
  0.5–2, ≥ 5 runs lie close, the held-out RMSE beats the prior-mean run and the T1 calibration (also on
  the 18 glaciers 5a never saw), and one real GloGEM run at the adopted parameters confirms it.
  Otherwise nothing is written back to IDL.
- Afterwards: per-glacier corrections δ for prediction and write-back; your old-versus-new trend test
  (Colle Gnifetti, Col du Dôme; dropping the Grenz pair, which compares two places, needs your OK); the
  temperate/cold classification.

## Not scheduled (documented)

| item | revisit |
|---|---|
| melt-season heat on ice (~+1.8 °C on the Hintereisferner tongue) — new physics | after the calibration works, only if ~2 °C remains, with your go-ahead |
| snow timing, leftover snow, depth reference, snow density | after Step 4; the depth reference in 5a |
| firn store; Herron–Langway stage 1; Tête Rousse's deep warm bias | together, at the next physics review |
| climate before 1942 (Colle Gnifetti); crevasse refreezing; strain heating; melting surface | documented / off |
| snowline calibration; winter snow > 2× off WGMS | your OK / with the mass-balance developers |
| Hintereisferner drift-correction stall | measured in Step 2 (c) |

## Verification

| step | check | pass |
|---|---|---|
| 0 | sine wave; last-day and first-month readings; folds | phase within 1 day; today's value; no shared model column |
| 0 | alps_r7 re-scored with 0a alone | best 2.337 °C at run 13 |
| 1 | profiles vs archived design 13; firn line logs | byte-identical |
| 2 | IDL test column: exact steady answer; Python copy | 0.0001 °C; 0.001 °C |
| 2 | firn bands; Hintereisferner / Grenz stop checks; drift vs 3000 cycles | 0.01 °C; within bounds; 0.01 °C |
| 4 | ensemble complete; table and keep rule | 32/32 runs, 66/66 entities; shown to you |
| 5 | notebook reproduces the alps_r7 numbers; `target = entity` = today | exact |
| Then | adoption checks; one real GloGEM run | all pass, or not adopted |

## Housekeeping

- **Order:** Steps 0 and 1 in parallel → Step 2 → (Step 3 decision) → Step 4 → stop → 5a → 5b →
  campaign. Runs one after another on vierzack03 (they share the GloGEM tree, compiled at each launch,
  and the glacier list in the catchments folder); no GloGEM edits while any run is going; IDL licence use
  checked before each launch; no process is ever killed.
- **Commits:** small, one per item, pushed once verified; reasoning in the message; book commits
  separate; never config files. GloGEM: Step 1 logbook, Step 2 fix. icetemp: Step 0 items, readers, test
  column, Step 4 scoring, notebook 06, Step 5b items.
- **Developers' model-changes document** (owned by the GMIP4 session; we hand over text, do not edit):
  Step 1 diagnostic files; Step 2 snow conduction fix (what was wrong, effect from Step 4; no effect on
  GMIP4, which runs with the thermal module off).
- **Book:** chapter 6 corrected (−4.4 to −5.0 °C; the cold tongue explained by the snow bug; firn memory
  and snowline decided); new chapter 7 "The snow, the comparison and the calibration target" with
  figures; this plan archived as `development_history/plans/2026-09-30_snow_comparison_target_plan.md`.
- **Memory:** firn-extent plan, f_ref/f_ins, snow insulation, comparison timing, target decision.
- **Storage:** raw profiles to `net_scratch/icetemp_profiles/`; home at 6.5 of 10 GB; the 456 MB of
  extracted alps_r7 profiles in the scratch plan folder deleted after 5a.

## Decisions for you

**Now**

| topic | recommendation | why |
|---|---|---|
| firn memory (Step 3) and snowline calibration (3c) | leave both out; record the findings in chapter 6 | their premises came from the wrong mass balance; no borehole or calibration sees the store |
| borehole l9b (Grenz tongue) | set to 2545 m as a documented override (or exclude it) | its position is published, only its elevation was estimated, and 2015 m lies below the glacier |
| where the snow fix acts | ice bands only | keeps your "snow inert on firn" |

**Later, when the step comes up**

| topic | when | recommendation |
|---|---|---|
| keep-only-what-pays thresholds | Step 4 | as in the Step-4 table |
| per-glacier offset δ in the fit | 5a | no; estimate it afterwards (with it held-out RMSE 2.75–3.12 °C, without 2.34–2.45) |
| error model | 5a | the August one (3.5 °C², 15 m) if it passes the misfit-ratio check; the fully fitted one as a sensitivity |
| RMSE margin for separation | 5a | 0.10 °C |
| if no target passes | 5a | two parameters with the measured cost; f_ref alone only if nothing separates |
| T2 or T3; readings at the melting point; Colle Gnifetti below 40 m; freeing f_adv | 5a | the pre-registered rule; shown both ways; out of the fit, reported apart; only if it passes the rule |
| campaign: f_ref range, scope | Then | 0.1–1.0; all 116 entities |
| Grenz tongue trend pair | Then | drop it (two places, not two times) |
| Hintereisferner's remaining ~2 °C | after Step 4 | document it; no new physics before the calibration works |
