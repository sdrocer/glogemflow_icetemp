# Snow insulation, a fair comparison, and a calibration target that tells two parameters apart

*(Revised 2026-09-30. Replaces the 2026-09-29 draft in this file; the firn-extent plan of
2026-09-28 is archived in `development_history/plans/2026-09-28_firn_extent_plan.md`. Every number
below comes from an investigation that an independent review re-checked; a last consistency check
re-read the run outputs and corrected a few. Where they disagreed, the table "Where the checks
disagreed" says which one this plan follows.)*

## Context — your three questions

### 1. "This is just after the spin-up, not after we calibrate?"

**Yes.** Every number in this plan comes from **forward runs**, each with one fixed set of
parameters, or from a Python copy of one model column (Step 2) or an approximate simulation
(Step 3). **No calibration has been run yet.**

```
 one forward run (about 11 min for the 7 test glaciers)
 ┌────────────────────────┐   ┌───────────────────────────────┐   ┌─────────────────────────┐
 │ spin-up to equilibrium │ → │ main run, Oct 1971 – Sep 2025 │ → │ compare with every      │
 │ in the 1942–71 climate │   │ same f_ref, f_ins, f_adv      │   │ borehole reading, at    │
 │                        │   │                               │   │ its own month           │
 └────────────────────────┘   └───────────────────────────────┘   └─────────────────────────┘
 alps_r7 = 32 such runs, with the parameters spread over their ranges (not fitted)
 "best run" = the lowest RMSE of the 32: design 13, f_ref 0.52, f_ins 1.09, f_adv 1.05

 calibration (after Step 5, not run yet): an emulator trained on such runs lets a sampler
 search *between* them for the best parameters and their uncertainty
```

A calibration can only move the parameters where they act. On ice sites all three move the
temperature today by just 0.1–0.3 °C per prior standard deviation, so no calibration can close the
4–5 °C Hintereisferner misfit. (Step 2 gives f_ins a hold on the ice through the snow, but about
2 °C of the misfit remains even then.) That is why a misfit is first fixed if it is a bug, and
documented if it is missing physics.

### 2. "Shouldn't we first check whether the recent changes improved the model? Focus on the calibration."

**Agreed. The plan is built that way:**
- **No new physics.** Only two kinds of change: bug fixes (the code, the comparison or the data do
  not do what they should) and measurements.
- **Firn memory is not in the model, and we recommend not adding it now.** We checked its premise
  first: with the correct mass balance it would change almost nothing the boreholes or the
  calibration see (Step 3). The same holds for the snowline calibration (3c of the 28 Sep plan).
  You had asked for both, so this is your call.
- **Step 4 measures** what the fixes bring, before anything else, on the same 32 parameter sets and
  the same readings. Anything that is not a bug fix stays only if it pays.
- **Step 5, the calibration target, is the main work.**

What the recent changes brought, each measured when it landed:

| change (2026) | kind | measured effect | now |
|---|---|---|---|
| conduction on the stretched grid, implicit solver, bed node, column down to the borehole (22–24 Sep) | bug fixes | old conduction 20× too steep in the top 10 m; 238 of 776 readings had been invisible | in |
| thermal spin-up: 30-year climate, main-run physics, true equilibrium (24–28 Sep) | bug fixes, plus model choices you approved (no meltwater heat in ice, lapse-rate rule, burial, relative depth) | RMSE 3.75 → 2.51 °C on 482 readings (Tête Rousse, Grenz/Gorner)¹ ² | in |
| each reading compared at its own month (28 Sep) | comparison fix | RMSE 3.22 → 2.50 °C, same run¹ | in |
| readings between two nodes interpolated (28 Sep) | comparison fix | 702 → 783 readings used; RMSE unchanged | in |
| each glacier's own (RGI7) mass-balance parameters (29 Sep) | input fix | ensemble best 2.37 → 2.34 °C, median 2.92 → 2.66 °C; ELA follows WGMS. At the old thermal parameters the fit first got worse (2.58 → 2.88 °C): they had been tuned to the wrong mass balance | in |
| strain heating, 900× unit error fixed (29 Sep) | process | fit worse at every depth | fixed, **off** |
| surface held at 0 °C in melt months (29 Sep) | process | RMSE 2.88 → 3.43 °C | **off** |

¹ measured with the wrong (RGI6-numbered) mass balance: valid as a before/after comparison only.
² scored on the model's final year, before each reading was compared at its own month.

Every bug fix stayed; both processes we tested were switched off because they did not pay. These
effects were measured one at a time, with different yardsticks, so they cannot be added up into one
before/after number. Step 4 is the first like-for-like test of a change on all 32 parameter sets.

**How much physics does a flowline model need?** Enough to get each band's thermal regime right,
because that sets the 10–50 m temperatures the boreholes and the calibration see. Beyond that, each
process brings a parameter that about 800 borehole points on 7 glaciers cannot pin down.

| process | why it matters for the boreholes | status |
|---|---|---|
| where the firn is | refreezing meltwater warms firn towards 0 °C; ice gets none | follows the modelled ELA, which matches WGMS on Hintereisferner since the RGI7 fix; the spin-up's starting firn area and Tête Rousse's firn history remain doubtful (Step 3) |
| surface temperature and its seasonal cycle | drives everything above ~15 m | from the mass-balance model's air temperature, T_off included |
| winter snow on the ice | the only insulation ablation-zone ice has | short-circuited → Step 2 (a bug fix) |
| ice flow: advection, burial, emergence | carries cold or temperate ice into the tongues | in, also in the spin-up |
| geothermal flux; a spin-up at equilibrium | set the deep ice | in |
| strain heating; surface at 0 °C in melt months | tested | off: both worsen the fit |
| meltwater heat in ice, crevasse refreezing, firn memory, climate before 1942 | local; each needs a new parameter or data we lack, or (firn memory) changes nothing the boreholes see | not now, documented |

The calibration can then correct what its two parameters reach: how much meltwater heat the firn
gets (f_ref), and how easily heat moves through firn and snow (f_ins). It cannot correct a process
that is missing at a few sites; those become documented limitations. A new process is worth
considering only when a misfit (1) is larger than the parameters can reach, (2) appears at several
glaciers, (3) has a physical cause with its own signature, and (4) pays in a paired ensemble test
like Step 4.

### 3. "Fit two parameters to reduce the borehole RMSE; they must be told apart for the chosen target."

**Exactly, and the target decides it** (Step 5):
- **Today's target** keeps one mean misfit per borehole entity. With it f_ins cannot be calibrated:
  the data alone leave it less certain than the prior already is.
- **Comparing per depth bin (T2), or per depth bin and month (T3),** tells the two apart. Estimated
  from the 32 runs, their correlation after a fit would be 0.57 or 0.34, instead of 0.87.
- **The catch:** today the lowest RMSE and the clearest separation put the parameters in slightly
  different places. The disagreement sits on Monte Rosa's high firn, and on glaciers left out of
  the fit it costs 0.05–0.15 °C of RMSE.
- **The rule for choosing** is fixed before we look, and applied to the Step-4 ensemble. In August
  you chose the whole profile shape as the target; T2 and T3 are that, and T1 stays as the reference.

**Words used below.**
- *Borehole entity:* one glacier's boreholes in one 10 m model band, e.g. Grenzgletscher@2505.
- *Point:* an entity's readings pooled into one 1 m depth bin; 769 on the 7 glaciers have a model
  value in all 32 runs.
- *Prior sd:* the spread allowed before any fit: f_ref 0.29, f_ins 0.33, f_adv 0.24.
- *Posterior:* the parameter values a calibration finds, with their uncertainty.
- *Held-out:* scored on a glacier system that the fit did not see.
- *Fold:* one round of "leave one glacier system out": fit on the others, score on it.
- *Emulator:* a fast statistical stand-in for GloGEM, trained on forward runs.

---

## What we found (what the plan rests on)

1. **The correct mass balance put the firn line in the right place, and the best fit improved
   slightly.** Icetemp `65de865` gives each glacier its own RGI7 parameters. The modelled ELAs
   now match the calibration, and Hintereisferner's follows WGMS (offset −68 m, r 0.83, slope
   0.97). Ensemble `alps_r7`:
   - best RMSE 2.36 °C on all 769 points (2.34 °C on the 689 shared with the previous ensemble, the
     number in chapter 6);
   - median 2.67 °C;
   - best runs at f_ref 0.38–0.54, inside the range.
2. **Two misfits remain at the best run** (design 13):
   - **Hintereisferner tongue ice (2555–2955 m): −4.4 to −5.0 °C**, 41 points, mean −4.7 °C. The
     draft's "−3.5 °C" also averaged in the two firn boreholes (3065 m −0.1 °C, 3155 m −0.5 °C);
     chapter 6 repeats it and will be corrected. With the Step-0 timing fix the tongue looks
     0.1–0.3 °C colder still (−4.7 to −5.1 °C, mean about −4.9 °C): its summer readings are then
     compared with mid-month values.
   - **Colle Gnifetti below 40 m: +1.3 to +1.7 °C**, the missing cold of the climate before 1942.
     Accepted and documented, as you decided.
3. **Snow on ice barely insulates, because of a solver bug.** The snow is stored half a layer too
   high, and every layer boundary gets the average conductivity of its two nodes. Thin winter snow
   (0.1–0.3 m w.e.) passes heat 5–7× too easily. → Steps 1–2.
4. **The comparison has two bugs, the held-out test leaks, and one borehole is almost certainly
   listed at the wrong elevation.** → Step 0.
   - Each reading is compared with the model at the *end* of its month: 1 °C RMS error in the top
     5 m.
   - Profiles dated only to a year are compared with January.
   - "Leave one glacier out" kept parts of the left-out glacier in training.
5. **The reason for firn memory is gone.** With the correct mass balance the Grenz tongue is never
   inside the firn area, and its false cooling has disappeared. → Step 3.
6. **Whether f_ref and f_ins can be told apart depends on the target.** Today the level of the
   profiles and their shape ask for somewhat different parameters. → Step 5.
7. **Also found, deferred (your decision).** `firnice_perm_depth.pro` multiplies the Herron–Langway
   stage-1 depth by the accumulation rate, which the 1980 paper does not do.

### Where the checks disagreed, and what this plan follows

In every row a review re-ran the numbers, and the plan follows it. Rows marked * come from the last
consistency check (30 Sep), which re-read the run outputs.

| topic | investigation or draft | review | why the review |
|---|---|---|---|
| Hintereisferner tongue misfit | −3.5 °C | −4.4 to −5.0 °C (41 points) | −3.5 included two firn sites |
| how much too easily heat crosses the snow | 8× (draft), 9× | 5–7× at 0.1–0.3 m w.e. | one consistent snow depth for both numbers |
| how much of the Hintereisferner misfit the snow fix removes | "most" | about half: +2.2 to +3.0 °C against −4.7 | matched to the readings' months and depths |
| where f_ins goes after the snow fix | "may drift above 1" | pulled **down** by the tongue, **up** by the firn | sign checked in the column copy |
| what f_ins does | "insulation acts in winter, in the top metres" | on firn bands it scales the whole firn column all year; on ice it acts only through the snow | read in the code |
| where the snow logbook is collected | the `glogem.pro` call site | inside the thermal model | already silent in the spin-up; no `glogem.pro` edit, fewer merge conflicts |
| test of the IDL edit | the logbook's two resistances agree | an IDL test column against an exact answer | the logbook's numbers agree by construction |
| design-13 values for the identity check | 0.5243 / 1.0890 / 1.0550 | 0.5245 / 1.0893 / 1.054977 | what the ensemble actually wrote |
| Tête Rousse bands at the reading dates | firn until 2021 | ice in 2010 and 2015; first flips in 1991 | firn line log re-read (verified) |
| spin-up check | Grenz/Gorner alone, ~75 min | 00777, 01225, 01230 + 03116, ~1.5 h | chapter-5 timings; Hintereisferner is the glacier we care about |
| firn store rule | empty it when the mass balance says ice | keep it while the band is on the glacier | emptying wipes decades-old firn after 1–2-year exits (95 cases) |
| firn store benefit at Tête Rousse | 3155 / 3165 m down to 0 … +0.5 °C | at best +1.2 °C at 3155 m | the deep warm bias is the firn history, not memory |
| spin-up reference period | 1961–90 | 1942–71 | run configuration (verified) |
| held-out cost of shape targets | +0.33 °C | +0.05 to +0.15 °C | the folds leaked (reproduced) |
| comparison timing | "no change needed" | end-of-month snapshot: 1 °C RMS in the top 5 m | read in the IDL writer (verified) |
| training campaign | 40 + 3 runs, some at f_ref < 0.1 | 40 + 20 in two stages, f_ref ≥ 0.1 | a 40 + 3 design fails the adoption check in 73–99 % of draws; runs below 0.1 break the emulator |
| error model | per-glacier offset δ plus a per-borehole offset | never stacked unfitted; δ kept out of the parameter fit | double counting; δ took the level and raised the held-out RMSE |
| * Tête Rousse boreholes the snow fix reaches | 3135 / 3155 / 3165 m (36 points) | every band except 3205 m (7 entities, 91 points) | the model had 3145–3195 m as ice in some years before their readings too |
| * reference for the Step-1 identity check | `temp_1m/10m/50m`, bedrock and firn line files of `alps_r7_w1` | the design-13 archive (profiles) and any worker's firn line log | `alps_r7_w1`'s folder holds the output of design 31, its last run; the archive holds only the profiles |
| * pooled RMSE expected after Step 2 | 2.16–2.22 °C | 2.16 °C if only Hintereisferner moves, up to about 2.3 °C if the other tongues and Tête Rousse warm too | recomputed from the design-13 points |
| * the remedy when "RMSE" fails | count Monte Rosa's high firn by its level only | on alps_r7 this costs more, not less (+0.18 °C held-out) | 2.52–2.54 against T1's 2.34–2.36 °C |

---

## The plan

```
 Step 0  compare model and boreholes fairly (Python, no IDL) ──┐ in parallel
 Step 1  snow logbook (IDL, no physics change)               ──┘
 Step 2  snow conduction bug fix + exact tests            one 11-min run, 1.5 h spin-up check
 Step 3  firn memory: a decision only (recommended: not built)
 Step 4  measure: the 32-run ensemble again, paired with alps_r7    ≈ 70 min on 6 IDL sessions
         ── stop: you see the table ──
 Step 5  calibration target: 5a analysis → rule applied → 5b code (Python)
 Then    two-stage training campaign (2 parameters) → KO calibration → adopt or not
```

Steps 0–4 use the 7 test glaciers (RGI7 00773 Bossons, 00777 Tête Rousse, 01112 Hohsaas, 01225
Gorner/Grenz, 01230 Monte Rosa east face, 02103 Corvatsch, 03116 Hintereisferner). Single runs use
the best ensemble point, design 13.

### Step 0 — compare model and boreholes fairly (Python only, no IDL)

*Why first:* every later number is measured with this ruler, and fixing it costs no IDL time.

**0a. Read the model at the reading's date, not at the end of its month.**

*Analogy.* A thermometer read on 10 July is compared with the value the model printed on 31 July.
Near the surface, three weeks of summer warming lie between the two.

GloGEM writes each monthly profile after the month's last time step (in
`firnice_temperature_model.pro` the write comes after the time-step loop), and the matcher compares
every reading of that month with this end-of-month profile.

```
     June                      July                      August
 ─────┼─────────────────────────┼─────────────────────────┼───▶ time
      │      reading on 10 July ▲                         │
      │                         └── today: compared with the model on 31 July
      ● model on 30 June ─── interpolate ──── ● model on 31 July   (fix: value on 10 July)
```

  T_model(reading) = T_before + f · (T_after − T_before)

| symbol | meaning |
|---|---|
| T_before | model at the reading's depth at the end of the previous month |
| T_after | model at the reading's depth at the end of the reading's month |
| f | fraction of the month elapsed on the reading's day (0 on the 1st, 1 on the last day) |

- Readings fall mid-month on average (median f = 0.56).
- The correction is 1.0 °C RMS in the top 5 m (up to 4.9 °C), 0.3 °C at 5–10 m, 0.1 °C at 10–15 m,
  and at most 0.03 °C below 15 m.
- On alps_r7 the best RMSE goes from 2.36 to 2.34 °C (run 13 stays best); the median from 2.67 to
  2.64 °C.
- The Hintereisferner tongue looks 0.1–0.3 °C colder afterwards (−4.7 to −5.1 °C): its readings
  were taken in summer, when mid-month is colder than month end.
- Left in, the bug would pull a calibrated f_ins down by about 0.3–0.37: 3.5–5 times the
  calibration's own uncertainty, and exactly in the near-surface shape where f_ref and f_ins differ.

**0b. Profiles dated only to a year or a range.** Their month comes from the start date, so today
they are compared with January.
- On the 7 glaciers: 16 profiles, 168 of 1722 month rows, all deeper than 10.8 m. The error is small
  (median 0.02 °C, at most 0.95 °C).
- Fix: average the model over the range, and add the model's spread over the range to that
  reading's error. Keep any range that overlaps the run window.

**0c. "Leave one glacier out" must leave the whole glacier out.** Two leaks:

```
   RGI7 01225 (Gorner/Grenz)                   │  RGI7 01230 (Monte Rosa east face)
   Grenz@4445   Grenz@4475   Grenz@4485        │  Grenz@4455   Grenz@4465
   └───────── the Colle Gnifetti saddle: boreholes 2–33 m apart in elevation ──────────┘
   leaving out 01230 still trains on its neighbours across the outline
```

- The folds were RGI outlines, and the Colle Gnifetti saddle straddles two of them.
- The pipeline's own validation folds by glenglat name. That separates Gornergletscher@2525 from
  Grenzgletscher@2525, which read the same model column.
- Fix: one fold per glacier system (01225 + 01230 = Monte Rosa), both in the notebook and in
  `validation.Validator`. On alps_r7 the held-out RMSE cost of the shape targets drops from
  +0.33 °C to about +0.05 … +0.15 °C (Step 5).

**0d. One borehole almost certainly at the wrong elevation (your decision).** glenglat lists l9b
(ryser2013, id 292) at 2015 m.
- The elevation is "estimated" from map contours, while the published position is 65 m from l9a
  (id 176, 2545 m), drilled the day before.
- 2015 m lies below the model's lowest band (2175 m), so its 10 points (28–267 m deep) are compared
  with the wrong ice.
- Recommendation: set it to 2545 m in the icetemp data handling, as one documented override (and
  tell the glenglat maintainers). Alternative: exclude it.
- No new IDL run is needed: the model already writes the 2545 m band for Gorner@2545, and that
  column reaches 318 m.

**Then re-score alps_r7** from its archived profiles, without IDL. That is the baseline for Step 4
and the reference for every expectation below. The set of points is fixed here; Steps 2–5 use
exactly these points.

*Why it is safe.* Only the comparison changes. GloGEM and its output do not, so nothing goes into
the developers' model-changes document.

*What could go wrong.*
- Within a month the model does not move in a straight line: it steps every 10 days, and snowfall
  and melt come in pulses. A straight line between two month ends is still much closer than the
  month end alone.
- The expectations in Steps 1–2 were computed with today's comparison. The *changes* they predict
  hold; the baselines shift slightly (Hintereisferner by −0.1 to −0.3 °C).
- A reading in the run's first month (October 1971) has no month end before it; it keeps today's
  value.

*Tests.* A synthetic sine wave is recovered to within one day of phase. A reading on the last day
of a month gives today's value exactly, and so does a reading in the first month. T2 equals the
count-weighted mean of T3 to 1e-12. In every fold, no model column (glacier, band) appears in both
training and test.

*Commits:* one icetemp commit each for 0a, 0b and 0c, and for 0d if you agree.

### Step 1 — a logbook of the snow the thermal model sees (no physics change)

*Analogy.* Before calling a sleeping bag broken, measure how thick it is, put one thermometer
outside and one inside, and note how much heat gets through.

Per band and calendar month (means over Oct 1971 – Sep 2025) the model writes:
- the snow on the band, and its depth from its mass;
- the surface temperature, and the ice just under the snow (also at 1 m and 2 m);
- how hard it is for heat to cross the snow-holding top layers: as the solver does it, and as
  layered snow and ice should give.

```
 one ice band in January (Hintereisferner tongue, 2615 m)
   node 0 (surface) ●  T_surf = min(0 °C, air)                        logged: T_surf
                    │  snow 0–1 m (250 kg/m³)                          logged: sno, z_snow
   node 1 (1 m)     ●
                    │  snow 1–1.17 m, then ice
   node 2 (2 m)     ●  first node under the snow                       logged: T_ice
                    │  ice
   R = thickness ÷ conductivity, summed from the surface to the first node under the snow
```

| logged | meaning | unit |
|---|---|---|
| sno, z_snow | snow at the start of the month, and its depth from its mass | m w.e., m |
| T_surf, T_ice | surface temperature; ice at the first node under the snow, end of month | °C |
| R_used | resistance of the snow-holding top layers, as the solver actually uses them | m² K W⁻¹ |
| R_series | the same layers as snow and ice one after the other: what they should give | m² K W⁻¹ |

Three small files per glacier:
- `firnice_snowclim_<id>.dat`: band × month.
- `firnice_snowmay_<id>.dat`: the snow on 1 May of every year, to compare with the WGMS winter
  balance.
- `firnice_snowbh_<id>.dat`: the monthly snow at the borehole bands. Step 5 needs it to know whether
  a shallow summer reading sat under modelled snow.

**Where the code goes.** Only in the thermal module's own files, so `glogem.pro` is not touched and
there are fewer clashes when this branch meets the GMIP4 code.
- The values are collected inside `firnice_temperature_model.pro`, next to its existing output, which
  is already silent during the spin-up. R_used is taken from the conductivities the solver actually
  used, not recomputed.
- The arrays go in `prepare_output_firnicetemp.pro`; the files are written after the firn line log in
  `write_firnicetemp_validation.pro`.

**What we expect.** From a Python copy of today's column (same grid, formulas and snow
bookkeeping), Hintereisferner tongue at 2615 m, in January:
- 0.30 m w.e. of snow, 1.17 m deep, under a −14.7 °C surface;
- ice at 2 m at about −9.4 °C, only 5 °C warmer than the surface (about 10 °C after Step 2);
- R_used about 1.3 against R_series about 7.3 m² K W⁻¹: **heat passes 5–6× too easily.**

On 1 May the copy carries about 0.7 m w.e. WGMS has band winter balances for Hintereisferner only
for five years between 2010 and 2015 (not 2012): 0.60 m w.e. at 2600–2650 m (range 0.17–0.90).

**What the logbook decides.** Not whether Step 2 happens (a bug is fixed anyway), but what we expect
from it:
- **0.2–0.4 m w.e. on the tongue in December–February:** a large effect (+2 to +3 °C) is expected.
- **Clearly less, around 0.1 m w.e.:** a smaller effect; the copy with half the snow gives +1.2 °C
  at 10 m. We re-run the copy with the logged snow and use that as the expectation, and to scale
  Step 2's stop checks.
- **No winter snow on the tongue:** the fix changes nothing there; we drop the Hintereisferner
  expectation.
- **1 May snow more than 2× off WGMS, thinner or thicker:** a mass-balance question. It is written
  down (chapter 7, and for the mass-balance developers) but not changed here: the mass-balance
  calibration is the standard one shared with GMIP4, and this plan does not touch it. Thicker snow
  than WGMS would also make the Step-2 warming larger (see "overshoot" in Step 2).

**The run** (`snowlog_w0`, about 11 min on vierzack03): the 7 glaciers at design 13, with
`alps_r7_w1`'s training configuration and only the output folder changed. The borehole list stays
alps_r7's too, even if 0d has landed, so the profile files keep their names. The parameters are set
with the ensemble's own override writer, which wrote 0.5245 / 1.0893 / 1.054977; hand-typed values
(0.5243 / 1.0890 / 1.0550) would break the identity check.

*Why it is safe.* It only writes, and the check proves it: every `temp_ID*.dat` is byte-identical to
the archived design 13 (`/itet-stor/jabeer/net_scratch/icetemp_profiles/alps_r7_w1/design_0013.tar.gz`,
65 profiles on 7 glaciers), and the firn line logs equal those of any alps_r7 worker (they do not
depend on the thermal parameters). The other summary files (`temp_1m/10m/50m`, bedrock) cannot be
compared: the archive holds only the profiles, and `alps_r7_w1`'s folder holds design 31, its last
run. GloGEM HEAD (`e9c0cca`, clean tree) is the code alps_r7 ran with, so the check should hold.

*What could go wrong.*
- IDL shares one variable space across files. `g`, `m`, `n`, `A` and `jg` are taken; the new names
  were checked to be free.
- In an ensemble the files would be overwritten at every design point, so the archiver keeps them
  with the profiles.
- If the identity check fails, we stop and find out why before Step 2; a clash in the shared
  variable space is the first suspect. The fallback reference is a fresh design-13 run at the
  parent commit.

*Commits:* one GloGEM commit (the logbook) and one icetemp commit (the readers and the archiving).

### Step 2 — let heat cross the snow as the model's own equations say (a bug fix)

**When it goes in, and in what order.** It is a bug fix, so it does not depend on the fit.
1. Build the IDL test column (Tests, a) and run it on today's code. It should miss the exact answer
   by 1–2 °C, as the Python copy does. If it does not, we stop and find out why the copy differs.
2. Make the fix; the test column must now hit the exact answer.
3. Run the 7 glaciers (b), then the spin-up check (c).

The logbook of Step 1 only sets expectations.

*Analogy.* Picture a wall: 1 m of insulation in front of 1 m of concrete. Heat must pass both, their
resistances add up, and the insulation decides how much gets through. The solver instead averaged
the two conductivities, as if the insulation were half concrete, so heat passed almost as through
concrete.

**What is wrong.** Temperatures sit on nodes (0, 1, 2 … 9 m, then 14 m …). Each node stands for a
slice reaching halfway to its neighbours (chapter 4), and heat flows through the gap between two
nodes. Three errors combine:
1. **The snow is stored half a layer too high.** The snow filling 0–1 m goes onto node 0, the
   surface node, whose temperature is set by the air. Its heat capacity is lost the same way.
2. **Every gap gets the average conductivity of its two nodes.** That is fine inside smooth firn,
   but wrong where snow meets ice. A partly filled gap is wrong too: snow and ice densities are
   averaged first, then turned into a conductivity.
3. **The snow depth uses 300 kg/m³**, while the snow layers get 250, 300 and 360 kg/m³. The snow the
   column holds differs from the model's snowpack by −17 % at 0.3 m w.e. and +21 % at 1.5 m w.e.

```
          what fills the gap           the solver uses now                  after the fix
  node 0 ●  surface (air temperature)
         │  1 m snow, 250 kg/m³        (0.15 + 2.13) / 2 = 1.14 W/m/K       0.15 W/m/K
  node 1 ●
         │  ice                        2.13                                 2.13
  node 2 ●
  a half-filled gap (0.5 m of snow over 0.5 m of ice):
         blended to 575 kg/m³ → 0.78, averaged with ice → 1.46              1/(0.5/0.15 + 0.5/2.13) = 0.28
```

How much too easily heat crosses the snow, from the surface to the first node under it (f_ins = 1):

| snow (m w.e.) | 0.1 | 0.2 | 0.3 | 0.45 | 0.6 | 1.0 |
|---|---|---|---|---|---|---|
| depth (m) | 0.40 | 0.80 | 1.17 | 1.67 | 2.14 | 3.21 |
| too easily by | 4.8× | 7.2× | 5.8× | 3.7× | 1.8× | 1.4× |

The thin midwinter snowpack is hit hardest; the factor jumps wherever the snow crosses a node.

**The fix.** In a gap that holds snow, the resistances add:

  Δz / k_gap = s / k_snow + (Δz − s) / k_col

| symbol | meaning | unit |
|---|---|---|
| Δz | thickness of the gap between two nodes | m |
| s | snow inside that gap | m |
| k_snow | Calonne (2011) conductivity × f_ins, never above the ice value | W m⁻¹ K⁻¹ |
| k_col | the column below the snow: Yen (1981) ice at the local temperature, ≈ 2.1 | W m⁻¹ K⁻¹ |
| k_gap | what the solver uses for that gap | W m⁻¹ K⁻¹ |

This is the textbook rule for layered material (Patankar 1980). Three consistency repairs come
with it:
- **(a) Snow depth** follows from the snow's mass and the densities its layers get.
- **(b) Heat capacity:** a node whose slice holds snow gets the heat capacity of what the slice
  actually holds.
- **(c) Emergence:** the snow-filled top nodes no longer rise 3.6× faster than the ice beneath them;
  they rise with the ice (0.01 °C).

**What stays the same.**
- Firn bands: their top already is snow-density firn, and snow stays inert there (`d072e0f`).
- No meltwater heat on ice bands (`97dfbf2`).
- The spin-up's drift correction: it sums the solver's own coefficients, so it follows by itself.
- The density curve, the mass balance, and T_off in the surface temperature.

**Why it is a bug and not new physics.** A column in steady state has an exact answer (surface
−30 °C, fixed snow, 0.3 W m⁻² of heat from below): across each material the temperature rises by
heat flow × thickness ÷ conductivity. For 0.15–0.45 m w.e. of snow the Python copy of today's scheme
misses it by 1.1–2.2 °C; the fixed scheme by 2·10⁻⁹ °C. The snow conductivity was introduced in
`1aef7f6` precisely so that snow would insulate; the solver never carried that out.

**What the column copy shows.** Hintereisferner tongue, 2615 m, 1942–1971 climate, the model's own
snowpack, at equilibrium, no ice motion, f_ins = 1:

| | now | after the fix | fine-grid reference | no snow at all |
|---|---|---|---|---|
| annual mean, 10 m | −6.3 °C | −3.3 °C | −3.2 °C | −5.3 °C |
| annual mean, 20 m | −6.0 °C | −3.0 °C | −2.9 °C | −5.0 °C |
| January, ice at 1 m | −11.0 °C | −5.8 °C | −5.8 °C | |
| June, ice at 1 m | −5.2 °C | −2.3 °C | −2.1 °C | |
| effect of f_ins at 10 m, per 0.35 of f_ins (about one prior sd) | 0.04 °C | 0.52 °C | | |

- **Today the snow even cools the tongue, by 1 °C.** In winter the snow is thin and short-circuited,
  so it barely insulates. In June it is deep and less short-circuited, so it still shields the ice
  from the 0 °C surface.
- Nearly all of the change comes from the conductance (+2.7 °C). The snow depth adds +0.3 °C; heat
  capacity and emergence change it by at most 0.02 °C each.
- After the fix the 1 m grid is within 0.15 °C of the fine-grid reference at 5–20 m (0.2 °C at 1 m
  in June): the grid is not the problem.

The result holds when the assumptions change (warming at 10 m):

| assumption changed | warming |
|---|---|
| none (the table above) | +3.0 °C |
| ice rising (emergence, b = −2.05 m w.e./yr) | +2.7 °C |
| half the snow / 1.5× the snow | +1.2 / +2.7 °C |
| a milder made-up climate (air −3 ± 7.5 °C; 0.2 / 0.4 / 0.8 m w.e. by April) | +0.7 / +1.6 / +1.7 °C |
| snow settling through the season (250 → 480 kg/m³, as in the April pits) | 0.3 °C more |
| top snow layer at 300 instead of 250 kg/m³ | 0.7 °C less |

**How much of the Hintereisferner misfit this explains: about half.**
- At the readings' own months (June–October 1976 at 2555 and 2615 m, March–August 1976 at 2955 m,
  March–October 1972 at 2805 m) and depths (0.3–14.6 m) the copy gives +2.2 / +2.4 / +2.4 / +3.0 °C
  at 2555 / 2615 / 2805 / 2955 m, against a mean misfit of −4.7 °C (about −4.9 °C after Step 0).
  The 2955 m value used the 2805 m forcing and is an upper value.
- Near the surface in summer (1–2 m) the change is only +0.9 to +1.5 °C.
- The rest most likely comes from the melt season. For context only: wet snow held at 0 °C would add
  about +1.8 °C at 10 m. That is new physics and not part of this plan.

**Who else changes.** Every borehole on a band that the model had as ice, with winter snow, at some
time before its reading (misfits at design 13, model minus boreholes, today's comparison):

| site (RGI7) | points | now | expected |
|---|---|---|---|
| Hintereisferner tongue 2555–2955 m (03116) | 41 | −4.4 to −5.0 °C | +2 to +3 °C; less at 2955 m, which takes ice from the firn band above |
| Grenz tongue 2015–2605 m (01225)* | 52 | −0.4 to +1.0 °C | warmer: up to +1.7 °C near the surface, less at depth (94 % of the Grenz and Gorner tongue points are deeper than 20 m) |
| Gorner tongue 2525–2555 m (01225) | 46 | −1.3 to −0.2 °C | warmer |
| Tête Rousse 3135–3195 m (00777) | 91 | +0.2 to +1.9 °C | warmer: 3135–3165 m were ice in most years since 2006, the others in a few years before their readings; may get worse |
| Hohsaas 3275–3345 m (01112) | 21 | −1.5 to +0.5 °C | warmer where the band was ice (from the 2022 update on) |

\* includes l9b, which Step 0d may move to 2545 m.

These are 251 of the 769 points. The other 518 sit on bands that were firn in every month before
their readings and must not move (test b).

*Rough estimate of the pooled RMSE* (design 13, today's comparison): if only Hintereisferner moved
as the copy predicts, 2.36 → 2.16 °C. If the other tongues warm mainly near the surface and Tête
Rousse and Hohsaas by about +0.5 °C, about 2.2 °C; if the tongues warm by +1.7 °C at every depth and
Tête Rousse by +1 °C, about 2.3 °C. Step 4 measures the real number.

**Why this matters for Step 5.**
- f_ins now reaches the ice sites (0.04 → 0.52 °C per prior sd at 10 m in the copy), where f_ref does
  nothing: a fingerprint of its own.
- **But** the ice sites stay about 2 °C too cold, from the missing melt-season heat. A too-cold
  tongue asks for a *lower* f_ins (more insulation); the firn sites, +0.4 °C too warm at design 13,
  ask for a *higher* one.
- So f_ins could partly stand in for missing physics. Step 5 checks whether firn and ice points want
  the same f_ins, and fits with and without the Hintereisferner tongue.

**Tests.**
- **(a) The IDL code itself** (mandatory: this is the test that proves the edit). A one-band test
  column in icetemp (`scripts/icetemp_analysis/snow_conduction_harness.pro`, with a Python check)
  sets the inputs the model reads and includes the model unchanged.
  - *Steady state:* surface −30 °C, 0.3 W m⁻² from below, fixed snow of 0.15 / 0.30 / 0.45 m w.e.,
    both solver branches. The exact answer is computed from the snow geometry, independently of the
    code. After the fix it must match to 0.0001 °C; today's code should miss it by 1–2 °C.
  - *A real year:* the Hintereisferner 2615 m monthly forcing, compared with the Python copy. It must
    match to 0.001 °C before and after the fix. This also proves that the copy behind every
    expectation above is faithful.
- **(b) The 7 glaciers:** `snowfix_w0` against `snowlog_w0`, same parameters, Step-0 comparison.
  - the "who else changes" table, site by site;
  - bands that are firn in every month unchanged to 0.01 °C (the spin-up's jump tolerance: the
    number of jumps is set by the largest jump over all bands, so firn bands can shift that much);
  - all 769 points: RMSE and bias, by firn / ice / unknown, with and without Hintereisferner;
  - spin-up jumps per glacier, and no new "did not converge";
  - logbook: R_used equals R_series on ice bands (the edit is wired in); January T(2 m) − T_surf
    goes from about +5 to about +10 °C.
  - *Stop checks* (an edit error can push either way): if Hintereisferner 2555 or 2615 m warms by
    less than +1 °C or by more than about +3.5 °C, or the Grenz tongue by more than +2 °C (beyond
    the copy's upper bound), the edit or the real snowpack differs from the copy. We compare the
    logbook and the test column before going on. If Step 1 logged less snow than the copy, the
    thresholds scale with the copy re-run on the logged snow.
- **(c) Spin-up.** The drift correction against 3000 plain cycles, with the fix, over all profiles
  and months: Tête Rousse, Grenz/Gorner and Monte Rosa east face (as in chapter 5), plus
  Hintereisferner; about 1.5 h in one IDL session. Pass: within 0.01 °C for 01225 and 01230.
  Hintereisferner already falls back to 200 plain cycles; we report how far that is from the
  reference. Tête Rousse alone would test nothing: all 25 of its bands are firn in the spin-up, so
  the changed code never runs.
  - It is launched directly on vierzack03 from the Step-2 run's own configuration, with only the
    spin-up settings changed (plain cycles, 3000, no early stop). It runs longer than the runner's
    time limit, and the `prun.sh` helper of chapter 5 copies an older configuration (1991–2020).
  - It runs after (b), not beside it: the two use different glacier lists (see Housekeeping).
  - If it fails for 01225 or 01230, Step 4 waits until the drift correction is fixed, because every
    Step-4 run uses it.

*Why it is safe.* Only ice bands that carry snow change directly, and the ice they pass down the
glacier; firn bands keep today's numbers to within the spin-up's 0.01 °C tolerance (test b checks
it). GMIP4 runs with the thermal model off, so they cannot see it. The IDL test column proves the
edit against an exact answer before any glacier is run, and the rule itself is the one the model's
own snow conductivity (`1aef7f6`) already implied.

**What could go wrong.**
- The Grenz tongue or Tête Rousse may end up too warm while Hintereisferner improves. The fix stays
  (a bug fix); the result tells Step 5 that f_ins is pulled two ways.
- **The fix could overshoot:** a tongue ends up too warm, or the pooled RMSE gets worse. It still
  stays. We first rule out an edit error (test column, logbook, stop checks), then compare the
  logged snow with WGMS: too much snow is a mass-balance question, written down and not changed
  here. What is left goes to Step 5. Too-warm tongues and the too-warm firn sites would then both
  ask for a *higher* f_ins, so the calibration can take part of it back; the rest is documented per
  borehole group, like Colle Gnifetti's deep ice.
- Now that snow really insulates, three simplifications start to matter. They are recorded, not
  changed, before Step 4:
  1. the thermal step keeps the start-of-month snow for the whole month; mid-month snow would warm
     the Hintereisferner sites by another 0.5–0.6 °C;
  2. last summer's leftover snow is removed only in January;
  3. the model measures depth from the snow surface, the boreholes from the ice surface. In
     June–July the copy still has 1.5–2.5 m of snow on the tongue, where the shallow 1976 readings
     were taken. Step 5 decides the depth reference, using the `snowbh` file.
- The drift correction may need more jumps: test (c) catches it.
- alps_r7, and every sensitivity number from it, is out of date after the fix; Step 4 reruns it.
- The explicit solver (not used in training) stays stable: its stability factor at node 1 goes from
  0.83 to 0.91 (0.97 at f_ins 2), still below 1.

*Commits:* one GloGEM commit (the fix) and one icetemp commit (the test column, the Python copy and
the exact test).

### Step 3 — firn memory: checked; recommendation: leave it out for now

*Forward runs: `firnrgi7_w0` (RGI7 mass balance; f_ref 0.87, f_ins 1.01, f_adv 0.94) against
`firnlog_w0` (the wrong RGI6-numbered parameters), plus the alps_r7 ensemble. The store's effects
come from an approximate simulation, because the model does not write per-band balances.*

**In short.**
- **The reason for this step is gone.** It was that the Grenz tongue cooled while its boreholes
  warmed. With the correct mass balance the tongue is never in the firn area, and the model warms
  there.
- **A store only matters where a band becomes firn for the first time.** In the whole run that
  happens at three updates only (1977, 1978 and 1981: 106 bands), as the firn area grows down from
  its spin-up extent (1942–71 climate) to the lower ELAs of the late 1970s. No borehole sits on those
  bands. The 234 later entries are bands that had been firn before: their old firn would still be
  there, so a store would leave them as today.
- **The calibration would not notice.** Tête Rousse and Hohsaas are the only glaciers with boreholes
  on bands that move in and out of the firn area. They hold 127 of the 769 points, but only 0.8 % of
  the information on f_ref and 0.2 % on f_ins.

**Recommendation: no firn store now.** Record the findings in chapter 6, and come back to the store
together with the meltwater depth (the Herron–Langway question), which decides how big a store would
be. This reverses your "do firn memory now", so it needs your OK. The snowline calibration of the
28 Sep plan (3c) goes with it: its premise, a firn line that reacts three times too strongly to
climate, came from the wrong mass balance, and with RGI7 the Hintereisferner ELA follows WGMS. If you
want the store anyway, a corrected minimal version is described below; it will not move the fit.

**1. The Grenz tongue.** Temperature change at fixed bands, 20–110 m deep, from 1979–80 to 2008–10:

| parameters | 2505 m | 2525 m | 2545 m | 2555 m | 2605 m |
|---|---|---|---|---|---|
| old, wrong (`firnlog_w0`) | −1.03 | −1.93 | −2.17 | −2.42 | −2.56 °C |
| RGI7 (`firnrgi7_w0`) | +0.55 | +0.52 | +0.46 | +0.42 | +0.36 °C |

- **Old parameters:** the firn area reached down to 2315 m, so the tongue started temperate and then
  cooled.
- **RGI7 parameters:** the firn area never comes below 2795 m. The tongue stays ice and warms.
- **The borehole pair is not a trend to match.** The 1979 borehole (2600 m, position estimated) and
  the 2008 boreholes (2553–2572 m) differ by +1.2 °C at 20–110 m, but also by +1.9 °C at 110–300 m,
  deeper than 29 years of heat conduction reach (about 35 m). The difference mostly reflects two
  different places.

**2. Where bands enter the firn area.** The thermal model starts from the firn area of its spin-up
climate, 1942–71. On five glaciers that firn area ends 115–172 m above the ELAs of the snowy 1970s,
and the 1977 and 1978 updates move its lower edge down by 170–260 m.

```
   Hintereisferner (03116), band 2805 m: warmest temperature of the year at 14 m depth
   1978   −6.9 °C   ice
   1979    0.0 °C   ← first year in the firn area: 49 m of firn, meltwater down to 30 m
   1982    0.0 °C
   1985   −1.8 °C   ← ice again after the 1982 update
   1990   −4.4 °C
```

| entries into the firn area, 1977–2025 | bands | when | warming at 14 m within two years (median) |
|---|---|---|---|
| for the first time | 106 | 1977 (50), 1978 (52), 1981 (4) | +3.1 °C (from −3.1 °C) |
| again (they were firn before) | 234 | 1979–2025 | +0.7 °C |
| neighbours that stay ice (comparison) | 546 | same years | 0.0 °C |

- **First entries are an artefact of the start.** In reality a band that turns firn carries a few
  metres of young firn on cold ice, and the meltwater refreezes in them.
- **Re-entries are real.** The old firn is still there and gets wet again, which is what the model
  does today.

**3. Emptying the store whenever the mass balance says ice would be wrong.** The firn area flickers
(2–4 bands per glacier switch each year). In 95 cases a band that had been firn for at least 10
years dropped out for only one or two years; firn that old holds tens of metres water equivalent and
cannot melt that fast. Emptying would put 298 band-years on thin new firn. Hohsaas 3275 m, for
example, was out in 2008 only: emptying would keep its meltwater in the top 1–6 m for 14 years, where
the model is already 1.2–2.1 °C too cold.

**4. Tête Rousse is too warm for another reason.** In the spin-up all 25 of its bands are firn, with
meltwater down to 32–67 m, and the firn line stays at or near the lowest band until 2000. The 2010
borehole at 3155 m reads −2.0 to −2.8 °C down to 46 m, which needs decades without meltwater. Even
the best case for any store (the column of the 3145 m band, which had not re-entered the firn area
by 2010) moves 3155 m only from +1.95 to +1.2 °C, and the readings at 26–46 m stay 1.4–1.9 °C too
warm. The misfit belongs with the meltwater depth (Herron–Langway) and Tête Rousse's balance
profile.

**If you want the store anyway: the corrected minimal version.**

*Analogy.* Each band keeps a savings account of left-over snow. Good years pay in and bad years
draw out, up to a ceiling equal to today's full firn column. The account stays open while the band
is on the glacier, even in years when the mass balance calls the band ice. Where the mass balance
says firn, the account sets how deep the meltwater can go.

```
   store S ▲  ceiling C = today's full firn column ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
  [m w.e.] │                                  ┌──────── full: exactly as today
           │                         ┌────────┘
           │        ┌────────────────┘  a bad year takes a little, not everything
           │   ┌────┘
           └───┴──────────────────────────────────────────────────▶ years
             first positive years: a few metres of young firn;
             the meltwater stays in them, and the ice below stays cold
```

| symbol | meaning | unit |
|---|---|---|
| b | the band's balance of the year, from the mass-balance model | m w.e. |
| S | firn store: the water equivalent of the band's firn | m w.e. |
| C | ceiling: firn mass from the surface down to the band's Herron–Langway depth (today's full column) | m w.e. |
| z_HL | Herron–Langway depth, today's meltwater depth | m |
| h | firn thickness: the depth at which the firn mass counted from the surface equals S | m |

```
  start of the run   S = C on bands that are firn in the spin-up; S = 0 on all others
  end of each year   every band on the glacier: S = min(C, max(0, S + b));  off the glacier: S = 0
  following year     where the mass balance says firn: meltwater down to the shallower of z_HL and h;
                     less than 1 m of firn (S below 0.28 m w.e.) counts as ice
```

| depth from the surface | 1 m | 5 m | 14 m | 30 m | 49 m |
|---|---|---|---|---|---|
| firn mass above it | 0.28 | 1.96 | 7.8 | 19.8 | 36.0 m w.e. |

- **Why it is safe.** Where firn is stays the mass-balance model's decision: one firn line. The mass
  balance, the firn line log, the spin-up and every band whose store is full stay as today. No new
  parameter, and no effect with the thermal model off (GMIP4).
- **What changes.** The 106 first entries keep cold ice under a few metres of young firn. A full band
  that loses mass in a bad year gets a meltwater limit about one layer shallower. Re-entries stay as
  today.
- **What could go wrong.** The heat of the 1977/78 entries no longer travels down the Gorner tongue,
  so its 2004–10 boreholes, already too cold, could get colder. Young firn follows the top of the
  density curve, so 1 m w.e. of it is 3 m thick instead of 1.7–2 m. A band returning after a long
  absence still gets almost its full column, because the Herron–Langway ceiling is 60–120 m on most
  firn bands.
- **Tests.** One 7-glacier run against the run just before it: the first month of every profile, the
  firn line log and the mass balance identical; S never above C; every reading within 0.01 °C except
  at Tête Rousse 3145 / 3155 / 3205 m and on the Gorner/Grenz tongue; the median 14 m warming of the
  106 first entries below 1 °C (today +3.1 °C). It stays unless a borehole group gets worse by more
  than 0.05 °C RMSE.
- **If it changes nothing** the boreholes can see (the expected outcome), that is itself the answer
  to your question 2: the extra code buys no fit. It then stays only if you want the start-up
  artefact gone for its own sake.
- **In Step 4** the store would run together with Step 2; its own effect is the single paired run
  above.
- **Cost.** About 15 lines in 6 GloGEM files, one 11-min run, one commit, one model-changes entry,
  one paragraph in chapter 6.

### Step 4 — measure what the changes bring, before anything else

*Analogy.* Weigh the same 32 samples on the old scale and on the repaired one. The parameter sets,
the borehole points and the comparison are the same; only the model changed. The difference for
each sample is the effect of the fix.

```
 design point        1    2    3   …   13   …   32      (the alps_r7 design, seed 11)
 alps_r7 (today)     ●    ●    ●        ●        ●   ─┐  both scored with the Step-0 comparison,
 alps_s2 (Step 2)    ●    ●    ●        ●        ●   ─┘  on the same points
                     └─ the difference at each point = what the fix does at that parameter set
```

**The runs.** `alps_s2_w0…w5` rerun the 32-point alps_r7 design with Step 2 in (and the store, if you
choose it) on the 7 glaciers, with raw profiles archived to net_scratch. About 70 minutes on 6 IDL
sessions on vierzack03 (alps_r7 took 67 min), after checking the IDL licence use. The GloGEM tree
must not change while they run: every run compiles it at launch. alps_r7 itself is re-scored from
its archive with the Step-0 comparison, without IDL.

**What we report** (one notebook section, and one table in the book):

| measure | why |
|---|---|
| pooled RMSE and bias on the same points: best run, and median of the 32 | the overall answer |
| change at each design point: median and range | is the gain the same everywhere in parameter space? |
| RMSE and bias per borehole group¹ | who gains and who loses |
| by firn / ice / unknown, and by depth (0–5, 5–15, 15–30, > 30 m) | where the change acts |
| strength of f_ref, f_ins and f_adv on firn points and on ice points, and how alike their patterns are | feeds Step 5a, and whether f_adv earns a place |
| where the best runs sit (f_ref, f_ins) | is the best point still inside the range? |

¹ Hintereisferner tongue ice and firn, Grenz tongue, Gorner tongue, Tête Rousse, Hohsaas, Monte Rosa
high firn, Colle Gnifetti below 40 m (the accepted limitation, shown apart), Mont Blanc (Col du
Dôme), Corvatsch.

**Keep only what pays.**

| kind of change | stays if | in this plan |
|---|---|---|
| bug fix: the code does not do what its own equations say, or the comparison is wrong | always. Its effect is reported. If the fit gets worse, something else was compensating (as the old surface offset did for the missing insulation, chapters 1–2): we look for it, we do not revert | Steps 0, 2 |
| artefact fix: removes a model artefact, adds no process | no borehole group gets worse by more than 0.05 °C RMSE | Step 3's store, if you choose it |
| new process or parameter | pooled RMSE better by at least 0.05 °C, both at the best run and in the median of the 32; no group worse by more than 0.1 °C; f_ref and f_ins not harder to tell apart | none scheduled |

**What we expect (estimates).**
- Hintereisferner tongue ice: from about −4.9 °C (after Step 0) to about −2.4 °C on average
  (−1.8 to −2.8 °C by site).
- Gorner tongue and Hohsaas better; Tête Rousse worse; the Grenz tongue warmer near the surface.
- Pooled RMSE on all 769 points: from 2.34 °C (after Step 0) to roughly 2.15–2.3 °C. It could also
  come out worse, if Tête Rousse and the Grenz tongue warm a lot (see "overshoot" in Step 2).
- f_ins on the Hintereisferner ice: from about 0.1 to about 0.5 °C per prior sd.

**Then a stop:** you see the table before Step 5's package code is written.

*What could go wrong.* The fit could get worse overall (Tête Rousse, Grenz). The fix stays; the table
then shows what was compensating, and Step 5 must know it. The three runs with f_ref < 0.1 stay in
for the pairing; Step 5 leaves them out of its slopes. If a run fails or loses entities (for example
a spin-up that no longer converges), the pairing uses the runs and points both ensembles have, and
the table says so.

### Step 5 — the calibration target (the main work)

*All Step 5 numbers come from the 32 alps_r7 runs (7 glaciers, 66 borehole entities, 769 points),
read through a smooth fit around the best run. No calibration has been run. Step 5a repeats
everything on the Step-4 ensemble before any package code is written.*

**The idea.** Think of the firn as a room that is heated in summer and cooled in winter.
- **f_ref is the heater:** how much meltwater heat goes into the firn when the water refreezes,
  filling the firn's cold content from the top down.
- **f_ins is the wall:** how easily heat moves through the firn and, after Step 2, through the
  winter snow on the tongues.

A warmer room can come from a stronger heater or a better wall, so the *mean* borehole temperature
cannot tell the two apart. Along this "level ridge" f_ref grows roughly with the square root of
f_ins, and between f_ins 1.0 and 1.9 the RMSE hardly changes (2.33–2.37 °C).

The two knobs leave different fingerprints in the profile:

```
 depth     refreezing (f_ref)                        conductivity (f_ins)
  0–5 m    summer heat, added from the top           how far the summer and winter swing reaches,
  5–20 m   a warm bulge down to where the water      and how late it arrives: less conductive firn
           runs out                                  damps it within a shallower depth, and delays it more
  > 20 m   only the long-term mean                   only the long-term mean
```

In the model, refreezing acts only on firn bands. There, f_ins scales the conductivity of the whole
firn column all year; the seasonal snow is not a separate layer on firn. On ice bands, f_ins acts
only through the winter snow.

**Three ways to compare the model with the boreholes (the "target").**

```
 one borehole entity (glacier, 10 m band)        what each target compares with the model
   depth   Jun-98  Nov-98  Jun-99
     2 m     x       x       x          T1 entity mean     1 number per entity          level
     5 m     x       x       x          T2 depth           1 number per 1 m depth bin   + shape with depth
    10 m     x       x       x          T3 depth + month   1 number per bin and month   + season
    20 m     x
```

T1 is today's target. T1 and T2 both average together readings from different seasons (26 % of the
depth bins) and from different years (25 %). Only T3 keeps them apart.

**What each target can tell apart** (estimated from the 32 runs, not from a calibration). f_adv is
fixed at 1.0. "Spread left" is the uncertainty expected after a fit divided by the uncertainty
before it: 1.0 means nothing was learned.

| | T1 entity mean (today) | T2 depth | T3 depth + month |
|---|---|---|---|
| numbers compared | 66 | 769 | 1674 |
| f_ref ↔ f_ins correlation expected after a fit | 0.87 | 0.57 | 0.34 |
| spread left, f_ref / f_ins | 0.26 / 0.72 | 0.11 / 0.28 | 0.09 / 0.22 |
| … with today's per-glacier offset δ in the fit | 0.46 / 0.91 | 0.13 / 0.29 | 0.11 / 0.25 |
| … without Monte Rosa (01225 + 01230) | 0.38 / 0.89 | 0.25 / 0.57 | 0.17 / 0.36 |

*With the Step-0 timing fix (re-run 2026-09-30): T2 0.60, spread 0.12 / 0.28; T3 0.35, 0.09 / 0.22.
The picture is the same.*

- **Today's target cannot calibrate f_ins.** The data alone would leave it wider than the prior
  (1.09); the 0.72 is the prior speaking.
- **Depth-resolved targets can.** That still holds with the emulator's own uncertainty included
  (T2 correlation 0.54–0.58, T3 0.34–0.37), at every candidate best point, and without the top 3 m.
- **Where the separation comes from:** the top 15 m of the firn sites (Monte Rosa above 3900 m,
  Bossons / Col du Dôme), not the ice sites (which leave 0.77–1.0 of the spread). Month-to-month
  changes carry a fingerprint of their own: f_ref and f_ins act on them almost independently
  (similarity +0.18), while on the depth pattern they act almost as one (−0.90).
- **f_adv** keeps 62–94 % of its spread in every target and is the weakest direction in T2 and T3.
  It stays fixed at 1.0 unless 5a, on the Step-4 ensemble, shows it has become useful (rule
  below).

**The catch: level and shape disagree, at a cost of 0.05–0.15 °C.**

```
 f_ins ▲
   1.9 │                        · · ○    ○ best f_ref for each f_ins (level only): RMSE 2.33–2.37 °C
       │                  · ○ ·
   1.1 │             · ●                 ● run 13, best of the 32: 2.34 °C (Step-0 comparison)
  0.85 │        ◆                        ◆ best T2 fit, Step-0 comparison: 2.44 °C
   0.6 │     ◇                           ◇ best T2 fit, today's comparison: 2.53 °C
       └──────────────────────────────▶ f_ref
           0.27 0.33     0.52      0.67
```

- **The disagreement is real.** Fitting the level alone and the shape alone gives optima 6.5
  standard deviations apart, before and after the Step-0 fixes. It sits on Monte Rosa's high firn,
  which holds 41 % of the points.
- **Its price is small.** On glacier systems left out of the fit (Step-0 comparison and folds) the
  RMSE is 2.34 °C with T1, 2.45 °C with T2 and 2.40 °C with T3, when choosing among the 32 runs; a
  smooth fit between the runs gives 2.34 / 2.48 / 2.39 °C. The earlier "2.44 against 2.77 °C" came
  from the fold leak. Six folds resolve such differences only roughly: the 90 % interval of T2 minus
  T1 is 0.02–0.16 °C (today's comparison).
- **An error model fitted to the misfits makes it costlier.** If the error model absorbs as much
  misfit as the residuals suggest (an offset per borehole, a 30–60 m depth correlation), it counts
  most of the level as model error; the shape then decides, and the held-out RMSE rises to
  2.56–2.60 °C. With today's per-glacier offset δ in the fit it is 2.75 °C (T3) to 3.12 °C (T2).

In plain words: for the lowest borehole RMSE the level must steer; to tell two parameters apart the
shape must steer. Today they disagree by a few tenths in f_ins. Steps 2–4 act on exactly this
near-surface shape and may shrink the disagreement, so we decide on the Step-4 ensemble.

**5a — analysis (a notebook; no package code, no IDL).**
`notebooks/06_calibration_target.ipynb` is built from this plan's prototype scripts and runs on the
Step-4 ensemble, re-matched with Step 0. For T1, T2 and T3 it reports the table above, the
level/shape split, the held-out RMSE one glacier system at a time (with a 90 % interval), and where
the information comes from (glacier, depth, season, firn or ice). It also checks:
- whether firn points and ice points ask for the same f_ins (Step 2's warning);
- the fit with and without the Hintereisferner tongue;
- readings at the melting point: 74 of the 769 points are within 0.1 °C of it, where the model cannot
  be warmer. Results with and without treating them as one-sided (your August question, still open);
- the depth reference for shallow summer readings (snow surface or ice surface), with the `snowbh`
  file;
- T3 without the top 3 m, because logger depths drift with burial or melt;
- Colle Gnifetti below 40 m (52 points, the accepted limitation), with and without; on alps_r7
  leaving them out changed nothing. Recommendation: keep them out of the fit and report them apart,
  so a known model limitation cannot steer the parameters;
- the three-parameter version of the table, to test your "f_adv fixed unless shown useful": after
  Step 2 the ice sites, where f_adv acts most, carry more signal;
- old against new readings (before and after 1979): if the old entities are clearly worse, your
  28 Sep option C (a model-error term for pre-1979 entities) comes back; otherwise nothing changes.

The decision rule, fixed before we look:

| check | passes if |
|---|---|
| told apart | f_ins spread left ≤ 0.5 and correlation ≤ 0.8, with two independent methods (a smooth response surface and the emulator), at the target's own best point and at 4 of the 5 best runs |
| consistent | fitting the level alone and the shape alone gives parameters within 3 standard deviations |
| RMSE | held-out RMSE ≤ T1's + 0.10 °C (the margin is yours) |
| honest error model | the misfit is about as large as the error model expects (ratio 0.5–2) |

Separately, f_adv is freed only if, with three parameters, its spread left is ≤ 0.5 and f_ref and
f_ins still pass "told apart"; otherwise it stays at 1.0.

What each outcome leads to:
1. **All four pass:** build T2, the smaller code change. Take T3 instead only if it cuts the f_ins
   spread by a further 20 % at no higher RMSE cost.
2. **Only "consistent" fails:** the parameters partly make up for a model error on Monte Rosa's high
   firn. Build the target if its RMSE cost is within your margin, and document the conflict like
   Colle Gnifetti's deep ice. Your call.
3. **"RMSE" fails:** you decide: two parameters with the measured RMSE cost, or f_ref alone with
   f_ins = 1.0 (which goes against your "at least two"). Counting the 19 Monte Rosa high-firn
   entities by their level only is tried, but on alps_r7 it does not help: it separates only just
   (correlation 0.78, spread 0.49) and costs more, not less (held-out 2.52–2.54 against T1's
   2.34–2.36 °C, +0.18 °C).
4. **Nothing tells them apart:** calibrate f_ref alone (against your "at least two"; your decision),
   and test parameters with a fingerprint of their own: a geothermal-flux scale, which acts only
   near the bed, or the percolation-depth factor after the Herron–Langway fix. A candidate first gets
   its own small ensemble test, like Step 4, and the campaign is re-planned for the parameters that
   remain.

On alps_r7 with the Step-0 fixes: T2 and T3 pass "told apart"; T3 passes "RMSE" (+0.05–0.06 °C) and
T2 fails it narrowly (+0.12–0.14 °C); both fail "consistent" (6.5); T2 passes "honest" (0.69), T3 not
yet (0.44: its seasonal error term must be re-measured). So alps_r7 would end in outcome 2 with T3,
provided T3's re-measured seasonal term passes "honest".

**5b — code** (only for the chosen target; icetemp; one commit per item). One switch,
`target = entity | depth | depth_month`; `entity` reproduces today's numbers exactly.

For T2:
1. *A likelihood with one misfit per depth bin*, with no averaging per entity:

   log L(θ) = −½ [ rᵀ Σ⁻¹ r + log det Σ ],   r = y − G(θ),   Σ = D + S + E(θ)

   | symbol | meaning |
   |---|---|
   | θ | the two parameters (f_ref, f_ins); f_adv = 1.0 |
   | y | measured temperature of each row: a 1 m depth bin (T2), or a bin and month (T3) |
   | G(θ) | the emulated model at the same depth and date |
   | r | the misfits |
   | D | measurement error: glenglat's value, ×3 for "estimated" profiles, plus the date spread of year-dated profiles |
   | S | model error: correlated with depth within a borehole (T2), plus a seasonal part between months (T3) |
   | E(θ) | emulator uncertainty (about 7–9 components, keeping 99.9 % of the variance) |
   | Σ⁻¹, log det Σ | weigh each misfit by the error expected there, and count correlated rows about once |

   S is chosen in 5a (see "Decisions for you") and then held fixed. The per-glacier corrections δ
   stay out of the parameter fit and are estimated afterwards (see "Then").
2. *Fixed parameters.* Today their handling works only when the fixed parameter is the last one
   (f_adv); a run with f_ins fixed would read the wrong bounds.
3. *Design over the free parameters only.* f_adv = 1.0 in every run (today's design code would still
   vary it), with a bound override for f_ref (0.1–1.0).
4. *Validation on GloGEM's emulator, not the analytical surrogate.* The calibrated parameters alone
   are scored through the emulator trained on GloGEM runs (today adoption is judged on the
   analytical surrogate), against the prior-mean run and the T1 calibration, one glacier system per
   fold. The adopted parameters are then confirmed with one real GloGEM run (see "Then").
5. *Configuration:* the target and the fixed parameters. The campaign yaml (yours to edit) needs
   `explained_variance: 0.999` (the field exists; its default is 0.99). At 99 %, half the f_ins
   information is thrown away (spread 0.46 instead of 0.26).

T3 adds: month rows kept next to the depth rows (so the baselines and the temperate/cold
classification still see one value per depth); per-(bin, month) model values with the Step-0
comparison; rows identified by position in the emulator, calibrator and validation (two months at one
depth must stay two rows); the seasonal error term from 5a.

*Tests.* `entity` gives today's likelihood bit for bit at 5 parameter sets. The fast and the full
likelihood agree. Two months at one depth give two emulator rows. In two-parameter mode f_adv never
moves. The chosen S gives a misfit ratio between 0.5 and 2 at the best run.

*Why it is safe.* 5a writes no package code, and 5b sits behind a switch whose default is today's
behaviour. Every decision is judged on glacier systems left out of the fit, with the fold noise
shown. GloGEM does not change, so there is no model-changes entry.

*What could go wrong.*
- **Seven glaciers, 32 runs.** Absolute numbers move with the method (T2's f_ins spread 0.21–0.53).
  The ranking T1 < T2 < T3 never changes, and the level/shape conflict shows in the real runs, not
  only in the smooth fit.
- **Monte Rosa dominates:** 41 % of the points, most of the separation, all of the conflict.
  Everything is reported per glacier and without it.
- **Thresholds.** Tête Rousse is at the melting point unless f_ref is almost 0. Runs below f_ref 0.1
  stay out of the slopes and out of the campaign.
- **The seasonal error term.** The first estimate had the wrong shape at depth: it let deep misfits
  lose their correlation within ten years, which the data do not. If the term is wrong, T3 over- or
  under-counts months; it is re-measured and checked with the misfit ratio.
- **f_ins has two jobs:** the firn's conductivity and the snow on ice. 5a checks whether both want
  the same value.

### Then — the training campaign and the calibration

*Analogy:* a survey. First a coarse grid over the whole area, then a fine grid where the target
turned out to be.

```
 f_ins ▲  stage 1: 40 runs over the whole box            stage 2: 20 runs around the stage-1 result
   2.0 │ ·    ·   ·     ·    ·                         │
       │   ·    ·    ·    ·     ·                      │           · ·
   1.0 │ ·   ·    ·   ·     ·    ·                     │         · ◉ ·
       │    ·   ·    ·    ·   ·                        │           · ·
   0.3 │ ·    ·     ·   ·    ·                         │
       └──────────────────────────▶ f_ref 0.1–1.0      └──────────────────────────▶ f_ref
```

- **Parameters:** f_ref 0.1–1.0 (prior and design), f_ins over its prior range (0.3–2.0), f_adv = 1.0.
  The three alps_r7 runs below f_ref 0.1 rank 27th, 28th and 32nd of 32 and break the emulator's
  error estimate.
- **Stage 1:** 40 runs spread over the box, an emulator keeping 99.9 % of the variance, and a first
  calibration with the chosen target. **Stage 2:** 20 runs within ±2.5 posterior standard deviations
  of the stage-1 result, a refit, and the final calibration.
- **Why two stages:** a single design of about 40 runs fails the existing check "at least 5 training
  runs close to the result" in 73–99 % of draws. With stage 2 it passes.
- **Scope:** all CentralEurope boreholes (116 entities on 25 glaciers), spin-up on, profiles
  archived, 6 IDL sessions on vierzack03.
- **Cost (estimate):** the 25 glaciers have about 2.6 times the bands of the 7 test glaciers, so
  about 25–35 min per run: roughly 3.5 h for stage 1 and 2 h for stage 2. One run is timed first:
  the pipeline stops a run after 30 min by default, so that limit is raised before the campaign.
- **Adopt only if:**
  - the sampling chains agree (R-hat ≤ 1.05), and the sampler's best point agrees with a brute-force
    grid search;
  - the misfit ratio is between 0.5 and 2;
  - at least 5 training runs lie close to the result;
  - the held-out RMSE (parameters alone, through the emulator, one glacier system per fold) beats
    both the prior-mean run and the T1 calibration; it is also reported separately for the 18
    glaciers that Step 5a did not see;
  - one real GloGEM run at the adopted parameters confirms the pooled RMSE. The emulator's own error
    (on alps_r7 about 0.3–0.4 °C per point) is as large as the differences being judged, so it
    cannot decide alone.
- **If a check fails:** nothing is written back to IDL, the current defaults stay, and the failing
  check says what to fix next.
- **Afterwards:**
  - per-glacier corrections δ, estimated from the remaining misfits, for prediction and the
    write-back to IDL;
  - the old-versus-new trend test you chose on 2026-09-28: Colle Gnifetti 1976–83 against 2003–21,
    and Col du Dôme 1973–80 against 2005–17. Dropping the third pair (Grenz tongue 1979–80 against
    2008–10) needs your OK: it compares two places, not two times. With T2 or T3 both epochs enter
    the fit, so the test is only partly independent, and it is reported as such;
  - the temperate/cold classification, reported as today;
  - the posterior correlation and spread, compared with 5a's prediction.

*What could go wrong.* A 25-glacier run can hit the time limit (hence the timed first run). The
18 new glaciers can bring misfits the 7 did not show; they are reported apart, and a new glacier
with a large misfit is looked at like Step 2 looked at Hintereisferner: bug first, then document.
Stage 2 can only sharpen what stage 1 found: if stage 1 lands at a prior edge, we report it and do
not adopt without your OK.

---

## Not scheduled (documented)

| item | what we know | revisit |
|---|---|---|
| melt-season heat on ice (wet snow at 0 °C, superimposed ice) | would add about +1.8 °C at 10 m on the Hintereisferner tongue (column copy); a new process | after the calibration scheme works, only if Hintereisferner stays ~2 °C too cold, and with your go-ahead |
| snow timing and depth reference | start-of-month snow held all month (+0.5–0.6 °C at Hintereisferner if mid-month); leftover snow removed only in January; depth from the snow or the ice surface | after Step 4; the depth reference is decided in 5a |
| snow density | firn curve 250–360 kg/m³ against April pits at 370–490 kg/m³; settling snow changes 10 m by only 0.3 °C, the midwinter top layer by 0.7 °C | with the item above |
| firn store (Step 3); firn outlasting the mass-balance firn area | no borehole or calibration effect now | with the Herron–Langway question |
| Herron–Langway stage 1 × accumulation (`firnice_perm_depth.pro`) | sets the meltwater depth, and any store's ceiling (60–120 m on firn bands) | your deferral; the next physics review |
| Tête Rousse deep warm bias (+1.4 to +1.9 °C at 26–46 m) | the modelled firn history and meltwater depth | with the Herron–Langway question and Tête Rousse's balance profile |
| climate before 1942 (Colle Gnifetti below 40 m, +1.3 to +1.7 °C) | documented limitation | — |
| crevasse refreezing (Gilbert et al. 2020) | backup; its signature is not seen | — |
| strain heating (unit fix `e9c0cca`), melting surface | both worsen the fit; off | — |
| snowline calibration (3c of the firn-extent plan) | premise gone: with RGI7 the ELA follows WGMS | only with your OK (Decisions, now) |
| winter snow amount, if Step 1 finds it more than 2× off WGMS | a mass-balance question; changes how much Step 2 can do | with the mass-balance developers; not in this plan |
| a model-error term for pre-1979 entities (option C of 28 Sep) | only if 5a shows the old entities clearly worse | 5a |
| Hintereisferner drift-correction stall (falls back to 200 plain cycles) | measured in Step 2 (c) | if the fallback is more than 0.05 °C from the reference |
| finer top grid (chapter 4, fix 2c) | after Step 2 the 1 m grid is within 0.15 °C of a fine grid at 5–20 m (0.2 °C at 1 m in June) | — |
| Grenz tongue 1979/2008 borehole pair | a difference between two places, not a trend | with the deep-ice questions |

## Verification per step

| step | check | pass |
|---|---|---|
| 0 | synthetic sine wave through the new comparison | phase within 1 day |
| 0 | reading on a month's last day, or in the run's first month | equals today's value |
| 0 | T2 from T3 | count-weighted mean, to 1e-12 |
| 0 | folds | no model column (glacier, band) in both training and test |
| 0 | alps_r7 re-scored with 0a alone | best 2.337 °C at run 13, as in the review's check; 0b and 0d then reported one by one |
| 1 | `temp_ID*.dat` against archived design 13 | byte-identical (65 profiles) |
| 1 | firn line logs against any alps_r7 worker's | identical (they do not depend on the thermal parameters) |
| 1 | logbook: tongue snow in Dec–Feb; 1 May 2010–15 against WGMS 0.60 m w.e. | read and noted (sets expectations) |
| 2 | IDL test column, steady state, both solvers | within 0.0001 °C (today 1–2 °C off) |
| 2 | IDL test column against the Python copy, Hintereisferner year | within 0.001 °C, before and after |
| 2 | bands firn in every month | within 0.01 °C |
| 2 | Hintereisferner 2555 / 2615 m; Grenz tongue | warmer by +1 to about +3.5 °C; Grenz by at most +2 °C (stop checks, scaled to the logged snow) |
| 2 | logbook R_used = R_series on ice bands | to print precision |
| 2 | drift correction against 3000 cycles (01225, 01230) | within 0.01 °C; 03116 reported; if it fails, Step 4 waits |
| 3 (if built) | first month, firn line log, mass balance | identical; S ≤ C |
| 3 (if built) | first entries' median 14 m warming | below 1 °C (today +3.1 °C) |
| 4 | ensemble | 32 of 32 runs, 66 of 66 entities |
| 4 | metrics table, keep rule applied | shown to you before Step 5b |
| 5a | notebook pointed at alps_r7 | reproduces this plan's alps_r7 numbers |
| 5b | `target = entity` | today's likelihood, bit for bit, at 5 parameter sets |
| 5b | fast vs full likelihood; two months at one depth; f_adv in 2-parameter mode | agree to 1e-9; two rows; never moves |
| Then | adoption checks | all pass, or the result is not adopted |
| Then | one GloGEM run at the adopted parameters | confirms the emulator's pooled RMSE |

## Housekeeping

- **Order.** Steps 0 and 1 in parallel; Step 2 after both; Step 3 is a decision. Then Step 4, a stop
  for you, Step 5a, the decision, Step 5b, and the campaign.
- **Compute.** All IDL runs on vierzack03, as for alps_r7 (the GMIP4 watchers leave it free).

  | step | IDL work | time | sessions |
  |---|---|---|---|
  | 0 | none: the 32 archived runs are re-matched in Python | minutes | 0 |
  | 1 | one 7-glacier run | ~11 min | 1 |
  | 2 | test column (minutes), one 7-glacier run, the 3000-cycle check on 4 glaciers | ~11 min + ~1.5 h | 1 |
  | 3 | only if built: one 7-glacier run | ~11 min | 1 |
  | 4 | 32 runs | ~70 min | 6 |
  | 5a, 5b | none | — | 0 |
  | campaign | 40 + 20 runs on 25 glaciers | ~3.5 h + ~2 h (estimate) | 6 |

- **Commits.** Small, one per step or item, and pushed once verified; the reasoning in the commit
  message, short code comments; book commits separate from code commits; never config files
  (`config.pro` and `config/*.yaml` stay yours).
  - GloGEM (branch `GloGEMflow`, main tree): the Step 1 logbook, the Step 2 fix, and the Step 3
    store only if you choose it.
  - icetemp: Step 0 (0a, 0b, 0c, and 0d if you agree), the Step 1 readers, the Step 2 test column and
    Python copy, the Step 4 scoring, notebook 06, and the Step 5b items.
- **Developers' model-changes document** (`_wt_v3/docs/source/development/model-changes-2026.md`):
  the GMIP4 session owns it, so we hand it the text (through you) and do not edit it. Only GloGEM
  changes get an entry (Steps 0 and 5 are Python). Proposed G24 bullets, handed in the page's own
  form (what was wrong, how found, effect, fix, test, submissions: none, since the module is off in
  GMIP4; main: has the older module, in which snow does not insulate at all):
  - Step 1: "new diagnostic files `firnice_snowclim/snowmay/snowbh` (firnice_temperature='y' only)";
  - Step 2: "seasonal snow on ice bands barely insulated after `1aef7f6`: it was stored half a layer
    high, and each layer boundary got the mean conductivity of its two nodes, so 0.1–0.3 m w.e. of
    snow passed heat 5–7× too easily; `<hash>` uses the series conductance of snow and ice, the snow
    depth from its mass and the snow's heat capacity; ice bands only; effect: Step 4's table";
  - Step 3: only if built;
  - two corrections to G24 as it stands: the thermal spin-up of the training runs uses the 30 years
    before the run (1942–71; 1961–90 is only the `settings.pro` default), and its count of commits
    that exist only on `GloGEMflow` grows by these.
- **Book.**
  - Chapter 6: correct "−3.5 °C" to −4.4 to −5.0 °C (41 points); add to its explanation of the cold
    tongue (a surface averaging about −5.5 °C) that the winter snow did not insulate (Step 2);
    settle the item "firn memory and snowline calibration to be re-decided"; add the Step 3
    findings.
  - A new chapter 7, "The snow, the comparison and the calibration target" (Steps 0, 2, 4 and 5),
    with three figures drawn by `figures/make_ch07.py`: the column copy (Step 2), the timing
    schematic (Step 0) and the level/shape ridge (Step 5); a row for it in the README's table of
    contents.
  - Once approved, this plan is archived as
    `development_history/plans/2026-09-30_snow_comparison_target_plan.md`.
- **Memory.** Update the notes on the firn-extent plan (firn memory decided), on f_ref/f_ins (the
  target analysis) and on snow insulation (the short-circuit); add notes on the comparison-timing
  bug and on the target decision; update the MEMORY.md index.
- **Storage.** Raw profiles go to `/itet-stor/jabeer/net_scratch/icetemp_profiles/<tag>/`; home is
  at 6.5 of 10 GB (measured 30 Sep). A run folder in `data/bayescal` takes about 18 MB for the 7
  glaciers. The 456 MB of extracted alps_r7 profiles in the scratch plan folder can be deleted
  after 5a.
- **Coordination.**
  - The GloGEM edits live on `GloGEMflow` in the main tree (clean, at `e9c0cca`). GMIP4 runs from
    the `_wt_v3*` worktrees with `firnice_temperature='n'`, so these edits do not change its
    results. We never touch those worktrees or the model-changes page.
  - IDL licences: 25 ETH-wide. GMIP4's last production campaign finished on 30 Sep, but more runs
    may follow (for example the Antarctic hindcast test), so we check the licence use before each
    launch.
  - Our own runs share two things, so they must not overlap: the GloGEM main tree (every run
    compiles it at launch, so no edits while any run is going) and the glacier list
    (`RGI11_CentralEurope_bayescal.dat` in the shared catchments folder, rewritten by each run's
    set-up). The 7-glacier runs, the 4-glacier spin-up check and the 25-glacier campaign therefore
    run one after the other.
  - `@`-included files are compiled even inside a false `if`, so a syntax error in the thermal
    module would stop every run compiled from this tree: begin/end balance is checked after each
    edit, and the test column compiles the file before any glacier run.
  - This plan kills no process. If one ever had to be stopped, liveness is re-checked right before,
    in a separate step.

## Decisions for you

**Now, before work starts**

| topic | recommendation | why |
|---|---|---|
| firn memory (Step 3) and snowline calibration (3c of 28 Sep) | leave both out for now; record the findings in chapter 6 | their premises came from the wrong mass balance; no borehole sees the store, and the calibration would not notice (under 1 % of the information); with RGI7 the ELA follows WGMS |
| borehole l9b, Grenz tongue | set its elevation to 2545 m, as a documented override; alternative: exclude it | its position is published, only its elevation was estimated; 2015 m is below the glacier's lowest band |
| where the snow fix acts | ice bands only; firn bands unchanged | the top of firn bands already has snow density; this keeps your earlier "snow inert on firn" |

**Later, at the step where it comes up**

| topic | when | recommendation |
|---|---|---|
| keep-only-what-pays thresholds | Step 4 | a new process must gain at least 0.05 °C of pooled RMSE and may cost no borehole group more than 0.1 °C; an artefact fix may cost no group more than 0.05 °C; bug fixes always stay |
| per-glacier offset δ in the parameter fit | 5a | no: it takes over the level (held-out 2.75–3.12 °C with it, 2.34–2.45 °C without); estimate δ afterwards for prediction |
| how much misfit the error model absorbs | 5a | the August error model (3.5 °C², 15 m depth correlation) if it still passes the misfit-ratio check on Step 4 (T3 adds its seasonal term); the fully fitted one as a sensitivity run: it lets the shape steer, and on alps_r7 it cost about 0.2 °C of held-out RMSE against T1 |
| RMSE margin to give up for separation | 5a | 0.10 °C |
| if no target passes | 5a | two parameters with the measured RMSE cost; f_ref alone only if nothing separates them |
| T2 or T3 | 5a | the pre-registered rule decides; T3 only if it clearly separates better |
| readings at the melting point | 5a | decide with the results shown both ways |
| Colle Gnifetti below 40 m in the fit | 5a | leave the 52 points out of the fit and report them apart (on alps_r7 this changed nothing) |
| when f_adv is freed | 5a | only if, with three parameters, its spread left is ≤ 0.5 and f_ref and f_ins still pass "told apart" |
| f_ref range of the campaign | Then | 0.1–1.0 |
| campaign scope | Then | all 116 CentralEurope entities (25 glaciers) |
| old-versus-new trend pairs | Then | keep Colle Gnifetti and Col du Dôme; drop the Grenz tongue pair, which compares two places |
| Hintereisferner's remaining ~2 °C | after Step 4 | document it; no new physics before the calibration scheme works |
