# Next steps (29 Oct): the free-air mass balance in the ice-temperature model, and back to the 2 m air's skill

*Approved by Janosch on 29 Oct 2026 (decisions A (a) and B (a); the bar in two steps; the full Step 5 redo only once
bar 1 is met).*

## Context

- **Your decisions today:**
  - **A (a):** version 1.0 takes the grid choice, f_ref 0.85 at f_ins 1.017.
  - **B (a):** the map's uncertainty for glaciers without boreholes is one massif's own spread.
  - The redo on the free-air mass balance follows the same two rules.
- **Your two requests:**
  1. Run the ice-temperature model on the mass balance recalibrated with the free air. It decides how much meltwater
     there is to refreeze in the first place.
  2. How good are we already below 15 m, as R² or RMSE?
- **Your goal:** the free-air model must get at least as good as the model on the normal (2 m) air was.
  - The free air stays in any case (standing rule: input data are chosen by how well they match their own
    observations).
  - So if we miss the bar, we fix the model, not the air.
- **The free-air mass balance is being made now** by the free-air chat window (its plan, steps 3 and 4).
  - It is recalibrating all 4,097 Alpine glaciers, our 23 included, on the free air. This started 11:15 on
    vierzack03 and should be ready about 13:30, in a new folder (`mbfreeair/step3/cal_fa`).
  - It also makes the same recalibration on the 2 m air (`cal_2m`), as a like-for-like reference.
  - The GloGEM code is the `freeair-mb` branch at 5ea744b. A last review of the thermal paths is running; that window
    tells me when the commit is final.
  - **Agreed split:**
    - That window does the recalibration and the check against GLAMOS (the Swiss measured mass balances).
    - I do the ice-temperature part (this plan).
    - Chapter 17 is mine. I hand that window the new f_ref and f_ins, and the borehole fit, for its chapter 15 and
      notebook.

## How good are we already? (your question 2)

Every number below is on the same **1,135 readings below 15 m**: each band at each depth, each hole and year, at 22
sites. "Unseen massifs" means each massif is predicted from the other six, as for a glacier without boreholes.

- **RMSE** is the typical miss, with big misses counting more.
- **R²** is the share of the spread in the measured temperatures that the model explains. The readings span 0 to
  −17 °C (spread 5.1 °C), so 1 means perfect and 0 means no better than one average value for all readings.

| model | how tested | RMSE | R² |
|---|---|---|---|
| **Step 5 (the plot in chapter 16)** | unseen massifs | **2.50 °C** | **0.76** |
| version 1.0 (grid choice 0.85 / 1.017) | all boreholes it was chosen on | 2.34 | 0.79 |
| version 1.0 (grid choice) | unseen massifs | 2.68 | 0.72 |
| **2 m air, three orientation parts** (o4, 23 Oct) | unseen massifs | **1.98** | **0.85** |
| 2 m air, five orientation parts (o5) | unseen massifs | 2.15 | 0.82 |
| free air, five parts, before f_sur (fb1) | unseen massifs | 2.71 | 0.72 |
| a straight line in elevation | unseen massifs | 3.67 | 0.48 |

Version 1.0 is on average 1.1 °C too warm, with a typical absolute miss of 1.6 °C. The best a band model can do
against single boreholes is about 1.7 °C (chapter 11).

**By height (unseen massifs):**

| | readings | measured spread | 2 m air (o4) RMSE / R² | free air now (grid) RMSE / R² |
|---|---|---|---|---|
| above 3800 m | 567 | 4.4 °C | 2.01 / 0.79 | **3.23 / 0.46** |
| zone 3000–3800 m | 390 | 1.8 °C | 2.20 / −0.50 | 2.21 / −0.51 |
| below 3000 m | 178 | 1.1 °C | (all near 0 °C; R² means little) | 1.33 |

**What this says:**
- **The R² of about 0.8 comes mostly from telling cold summits from temperate tongues.** The model does that well.
- **The free air costs 0.6–0.7 °C overall, and all of it above 3800 m** (2.0 → 3.2 °C there). This is the deep ice,
  30–200 m, at the high Mont Blanc sites (Taconnaz, the summit, Bossons) and Fiescher: too warm (chapters 13–14).
- **In the zone, R² is negative in every version, 2 m air included.** The model cannot yet order the zone bands from
  warm to cold. A straight line in elevation does better there (RMSE 1.37). That is a separate, older problem; the air
  did not cause it.
- **Five orientation parts instead of three cost 0.27 °C on the 2 m air** (1.86 → 2.13). On the free air the
  difference is only 0.08 (2.71 against 2.79).

## The thermal structure in general, by height (your comment)

**Two views beside the RMSE:**
1. **The regime class of each band:** cold, polythermal or temperate. The 105 bands are judged by the same rule as in
   chapter 16.
2. **Each profile split into two parts:**
   - its **level:** the mean temperature of one hole in one year;
   - its **shape:** how the temperature changes with depth around that mean.

**Bands in the right class:**

| | measured classes | v1.0 (0.85 / 1.017) | 2 m air (o4) | "everything is cold" |
|---|---|---|---|---|
| above 3800 m | 53 cold, 6 polythermal, 2 temperate | 53 of 61 | 53 of 61 | 53 of 61 |
| zone 3000–3800 m | 18 cold, 7 polythermal, 2 temperate | 10 of 27 | 11 of 27 | **18 of 27** |
| below 3000 m | 10 cold, 5 polythermal, 2 temperate | 9 of 17 | 10 of 17 | 10 of 17 |

(These counts are mine, approximately: one tolerance per band, 72 of 105 against chapter 16's 73.)

**Profile level and shape:**

| | | v1.0 | 2 m air (o4) |
|---|---|---|---|
| **above 3800 m** | level: miss / correlation | 2.67 °C / 0.83 | 2.22 / **0.89** |
| | shape: correlation | 0.63 | **0.84** |
| | change with depth (measured −1.27 °C per 100 m) | −0.37 | **−1.00** |
| **zone** | level: miss / correlation | 1.78 / 0.13 | 1.79 / 0.06 |
| | shape: correlation | −0.42 | −0.27 |
| | change with depth (measured +1.86 °C per 100 m) | +0.16 | +0.20 |
| **below 3000 m** | level: miss | 1.01 | 0.93 |

**What this says:**
- **Above 3800 m, the class is right in every version, but only because almost everything there is cold.** The model
  misses all 8 warm high bands (the Monte Rosa 1999 bands, Breithornplateau). So does "everything is cold".
- **The free air's loss above 3800 m is in the profile's shape.** Measured, the ice gets colder with depth (−1.27 °C
  per 100 m): the memory of a colder past.
  - The 2 m air reproduced that (−1.00); the free air makes the profile almost flat (−0.37).
  - This supports the earlier inference (chapters 13–14): the cold winter bias of the 2 m air stood in for the colder
    past that our spin-up lacks.
- **In the zone, no version predicts anything.**
  - Wrong level: the correlation is about 0.1.
  - Wrong shape: the correlation is negative. Measured, the zone warms with depth; the model is flat at 0 °C.
  - Wrong classes: worse than calling everything cold.
  - The air does not matter there. It is the same with the 2 m air.

## Why the zone is still so poor (my reflection)

**1. The model makes the zone all-or-nothing.**

At f_ref 0.85, 8 of the 9 zone site groups sit within 0.2 °C of 0 °C. At f_ref 0.03 several jump to −7 or −8 °C.
The measured values lie in between:

| zone site | measured | v1.0 (f_ref 0.85) | f_ref 0.03 |
|---|---|---|---|
| Taconnaz 3680–3830 m (to 113 m deep) | −6.7 °C | −0.2 | −7.9 |
| Taconnaz 3410 m | −2.6 | −0.0 | −7.1 |
| Alphubel South | −1.9 | −1.5 | −4.7 |
| Corvatsch | −1.8 | −0.0 | −3.3 |
| Tête Rousse | −0.9 | −0.1 | −0.2 |
| Jungfraufirn, Hohsaas, Altels | −0.2 to −0.7 | −0.0 to −0.2 | −0.1 to −8.3 |

**Why it behaves like a switch.** Refreezing meltwater carries a lot of heat. A rough estimate: 1 cm of refrozen
water warms 1 m of firn by about 3 °C. Even a small share of a summer's melt therefore heats the whole firn column to
0 °C.
- Our model puts each month's refreezing heat into the firn at once, down to the permeable depth.
- So the column flips between "soaked and temperate" and "dry and cold".
- Chapter 8 saw the same flip: a fixed 4 m meltwater limit made everything cold, and a limit that depended on the
  model's own state flipped back and forth.

The real zone has partly warmed firn. Water percolates unevenly (pipes and ice layers), and the winter cold
reaches different depths.

**2. Whether firn is there at all decides it, and the model's firn reaches too low.** The zone is the firn line
itself.
- Where the model has firn but reality has cold ice, the meltwater heats ice that, in reality, sheds it.
- Corvatsch is the proven case: with the observed firn and S, all 4 bands become right (chapter 11).
- *Inferred, not tested:* the mass balance before 2000 is 0.5–1 m w.e. per year too positive. That would keep too much
  firn in the spin-up and the early run.

**3. The air does not matter here, by construction.** Once the firn is soaked, it sits at 0 °C, however cold the
winter. That is why the 2 m and free air give the same zone.

What can move the zone is how much meltwater there is (the new mass balance's DDF and snow cover) and where firn is.
That is exactly what your step 1 changes, so the grid runs will show whether it helps.

**4. Local causes a band model cannot see.**
- The zone has only 7 sites, and they differ in kind: a shaded north face (Taconnaz), saddles near 0 °C (Jungfraufirn,
  Hohsaas), a small glacier whose drilled ice is 66 m thick where the model has 24 m (Tête Rousse).
- The measured spread in the zone (1.8 °C) is about the size of the 1.7 °C floor of a band model. So R² is the wrong
  measure there.
- Tête Rousse alone holds 319 of the 390 zone readings, so the per-reading numbers mostly describe Tête Rousse. The
  zone is better judged by classes and by site.

**What follows for the plan:**
- **The grid runs on the new mass balance (step 4) are read for the zone separately:** classes, and each zone site's
  level against the table above. Do the sites move off 0 °C?
- **After bar 1, a zone plan of its own** (for your approval). The candidates, in your order:
  1. **The firn extent:** the observed firn where it exists, and the snow-line product elsewhere.
  2. **The mass balance before 2000:** it is too positive, which may leave too much firn.
  3. **A gradual meltwater heating instead of the switch.** It would be tied to the climate, for example the winter's
     cold content against the summer melt, as chapter 8 required.

  The first two are parameters and inputs; the third is a change to a process, so it comes last.

## The bar (your choice: both, in two steps)

All scores use the score of chapter 16: each massif predicted from the other six, every outline counting equally,
with the grid choice of f_ref.

| | overall | above 3800 m | zone |
|---|---|---|---|
| **Bar 1:** the 2 m air with the same five parts (o5): a fair test of the air alone | 2.13 °C | 2.33 | 2.51 |
| **Bar 2:** the best 2 m air version (o4, three parts) | 1.86 | 1.87 | 2.52 |
| today: free air, five parts, f_sur 1.4 (option A) | 2.83 | 3.43 | 2.51 |

**"As good as" (my proposal):** at most 0.1 °C above the bar overall and above 3800 m, and the zone not worse by more
than 0.2 °C. These are the same tolerances as rule 9 and the earlier "no harm" rules.

## The steps

```
 today                                                                         tonight / tomorrow
 0. records ─► 1. set-up ─► 2. new mass balance ─► 3. new orientation ─► 4. 45 grid runs ─► 6. decision
 (A, B, the    (runner       arrives ~13:30;        factors (~1.5 h)      (~2 h): is the       bar 1 met: Step 5
  yardstick)    switch,      my checks                                    gap closed?          redo overnight
                gate test)                          5. what happened to the refreezing         not met: plan to
                                                       (alongside, no extra runs)              find the cause
```

### 0. Records and the missing check (today, about 1 hour)

- **The decision log** gets two rows: A (a) and B (a).
- **Chapter 16 and the rules plan's result section:** "decisions pending" becomes "decided".
- **The check that was still missing for B.** The ranges are tested around the grid choice instead of around Step 5's
  results. This uses the emulator only, takes minutes, and needs no GloGEM runs.
  - The centres are each fold's grid choice (1.0, 0.7, 0.85, 1.0, 0.55, 0.85, 0.85, at f_ins 1.017).
  - The spread is one massif's own: f_ref 6 × 0.16 ≈ 0.96. For f_ins it is Step 5's 1.2, since the grid holds no f_ins
    information.
  - It goes in a new script, `scripts/icetemp_analysis/step5/coverage_grid.py`. It imports the functions of
    `coverage_spread.py` (`draws`, `inside`) and writes `benchmark/coverage_grid.json`.
  - It is reported in chapter 16.
- **The comparison table above as a script**, so the numbers can be rebuilt:
  `scripts/icetemp_analysis/step5/yardstick.py`.
  - It reads `week3/freeair/score/points.csv` (o4), `week3/freeair/score_fb1/points.csv` (o5, fb1) and
    `benchmark/benchmark_points.csv` (option A, Step 5).
  - It writes `yardstick.json`.
  - Later runs add their own rows.
- **Commit and push** to glogemflow_icetemp.

### 1. Set-up for the new mass balance (while the recalibration runs, about 2 hours)

1. **The runner learns the mass balance's air.** A new setting `air_temperature` goes into the GloGEM config that the
   runner writes.
   - The default is `era5_2m`, which leaves today's config byte-identical.
   - The runner also copies `calibration/air_temperature.txt` along with the calibration files. Without it, GloGEM
     stops: that is its guard against using the wrong calibration.
   - Files:
     - `src/icetemp/calibration/runner.py`: the config template near line 245, the settings near line 406, the copy
       near lines 558–575;
     - `src/icetemp/calibration/config.py`;
     - `src/icetemp/calibration/pipeline.py:116`.
2. **A new config:** `config/step5_centraleurope_parts5_mbfa.yaml`, a copy of `step5_centraleurope_parts5.yaml`. Only
   these settings change:
   - the calibration folder (`cal_fa`);
   - `air_temperature: freeair`;
   - the new orientation-factor folder (step 3);
   - the GloGEM tree: `GloGEM_wt_mbfreeair` at the final commit;
   - the output: `data/bayescal/step5_parts5_mbfa`.

   Everything else stays as in Step 5: f_sur 1.4, five parts, FRHS, f_adv 1.0, the same readings and the same error
   model.
3. **New chain scripts**, copies of `chain_step5_bench.sh` / `chain_step5_points.sh`. They point at the new tree and
   commit and check them before every start. The old chains are not edited.
4. **Gate test:** one glacier, two runs.
   - With the switch off, the generated config is identical to today's.
   - With it on, GloGEM's calibration check passes and its log shows the free-air mass balance.
   - One five-part run of one glacier must finish OK.
5. **An independent reviewer** checks the runner change.
6. **The rules commit of the redo:** the runner change, the new config and the chains, committed and pushed. The rules
   are the same as in Step 5 (3467851); only the pointers change. The new chains check this commit.

### 2. The new mass balance arrives (about 13:30)

- I start nothing until the free-air chat window confirms three things: the thermal-paths review is done, the GloGEM
  commit is final, and `cal_fa` is checked.
- **My checks on `cal_fa` for the 23 calibration glaciers:**
  - all reached their target (standard tolerance, about 0.1 m w.e.);
  - `air_temperature.txt` says `freeair`;
  - a table of T_off, c_prec and DDF, old (Matthias' calibration, used by Step 5) against new.

### 3. Rebuild the orientation factors (about 1.5 hours on vierzack03 + 07)

The five-part factors were matched to each band's melt, and that melt changes with the new calibration. So they are
rebuilt.

- One melt-logging run of the 23 glaciers on `cal_fa` (FRHS, `--write-melth`, about 25 min).
- `part_factors.py` builds the factors into a **new** folder, `week3/orient_parts5_mbfa` (about 40 min).
  `orient_parts5_fa` stays untouched.
- **Check:** how much the factors at the borehole bands move against the old ones.

### 4. The 45 grid runs (about 2 hours on 03 + 07)

- **The same 45 runs as Step 5's benchmark**, on the new mass balance:
  - f_ref 0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85 and 1.00 at f_ins 1.017;
  - plus the starting values (f_ref 0.55, f_ins 1.0);
  - five parts each, under new run names and folders.
- These runs are also the benchmark runs of a later Step 5 redo, so they are not lost either way.
- **The parts are merged and scored with the same score.** Reported:
  - overall, the zone and above 3800 m, as RMSE and R²;
  - the 15–20 m and 30–200 m errors above 3800 m;
  - each massif;
  - the regime classes and the frozen beds;
  - the best f_ref per massif (does the Mont Blanc / Monte Rosa split move?).
- **The verdict against bar 1 and bar 2**, plus how far the best f_ref moved with the new mass balance. That tests your
  caution that f_ref carries the mass-balance calibration.

### 5. What the new mass balance did to the refreezing (alongside step 4, no extra runs)

**Per borehole band, old against new**, from the melt-logging runs (the old one exists from the free-air factor build)
and the grid runs:
- the melt that can refreeze (the thermal melt: free air without T_off, but with the new DDF);
- snowfall and burial (c_prec);
- the months on snow, firn and ice (the surface type comes from the mass balance);
- the firn extent;
- the refreezing heat.

This is the plain answer to "how much refreezing is there in the first place". The mass balance's T_off does not
enter the thermal melt. The new DDF, c_prec and snow cover do.

**Only if the score moves by more than 0.2 °C:** the same 45 runs on `cal_2m` (our own recalibration on the 2 m air,
about 2 hours). They tell apart what comes from the free air in the mass balance and what comes from our recalibration
differing from Matthias' (on a test glacier, T_off 0.00 against 0.12).

### 6. Decision point (with you, about 17:00–17:30)

- **Bar 1 met:** the full Step 5 redo runs tonight, with the rules commit of step 1.6. That is about 8–10 hours:
  - 40 design points and 20 more near the result, five parts each;
  - the seven left-out-massif fits;
  - the 40 runs at their results;
  - the 45 grid runs of step 4 as its benchmark.

  Then tomorrow: the results, rules 9 and 10 as fixed, chapter 17, and f_ref / f_ins old against new for the free-air
  chat window.
- **Bar 1 not met:** no Step 5 redo yet, since Step 5 must run on the final physics. Instead I write a plan for your
  approval to find the cause above 3800 m.
  - **The candidates already on record**, in your order (parameters within physical ranges first, new processes
    last):
    1. **The colder past in the spin-up.** Deep ice at 30–200 m remembers colder centuries that the 1942–71 spin-up
       climate lacks. The test at −2 K had the right shape, but acted only at 100–200 m: the 1971–2025 run resets the
       top 50 m (chapter 14).
    2. **Too much summer melt at Mont Blanc.** June–September melt is 1.5 times higher there, from its climate cell's
       spread of daily temperatures within each month. Test: the same spread for all glaciers (chapter 14).
    3. **Melt from an hourly energy balance instead of degree days**, at warm and cold bands. This is a new process,
       so it comes last.
  - **Already ruled out as the main cause:**
    - the surface offset f_sur: it shifts everything evenly and cannot separate Mont Blanc from Monte Rosa;
    - advection: at most 0.17 °C;
    - the mass balance's T_off (rank correlation 0.06);
    - snow burial: 0–0.8 °C.
  - Each test gets a reading fixed before its runs, as before.
- **Bar 2 (1.86 °C) comes after bar 1.** The free air moves the five-part / three-part difference to only 0.08 °C, so
  bar 2 hangs mostly on the same deep high-Mont Blanc ice.

## Coordination and cluster

- **The free-air chat window** owns the recalibration and the GLAMOS check. I pin no GloGEM commit before its message.
  Its 5 single-core past runs on vierzack03 after the calibration are off my path.
- **Before every launch on 03 or 07**, I send a message to that window. GMIP4 holds 04 and 06 until about 16:00.
- **Chains:** OMP_NUM_THREADS=1, nice 10, launches 45–60 s apart, launchers started on vierzack07. A running chain is
  never edited.
- **Untouched:** config.pro, `config/bayescal_centraleurope_25glaciers.yaml`, the deleted
  `notebooks/fig_03_CP_formula.pdf`, `firn_model/`. Chessjen stays excluded (01161).

## Records

- **Chapter 17:** the free-air mass balance in the ice-temperature model, with figures from `figures/make_ch17.py`.
  It covers the yardstick table, the grid result against the bars, the refreezing old against new, and then the redo
  or the cause plan.
- **The decision log and memory**, updated as results land.
- **Commit and push** each verified batch to glogemflow_icetemp.

## Verification

- **Runner:** switch off gives a byte-identical GloGEM config (dry run, compared). Switch on: GloGEM's calibration check
  passes and the log shows the free-air mass balance.
- **`cal_fa`:** 23 of 23 glaciers at their targets; `air_temperature.txt` = `freeair`.
- **Factors:** complete for every band and part; the borehole-band factors are compared with the old ones.
- **Grid runs:** 45 of 45 OK, every part merged, 1,135 readings scored.
  - Before scoring them, `yardstick.py` must reproduce today's numbers from the old point files: o4 1.855, o5 2.126,
    fb1 2.791, option A grid 2.831.
- **`coverage_grid.py`:** run with Step 5's centres and spreads, it must reproduce `coverage_spread.json`.

## Glossary

| term | meaning |
|---|---|
| RMSE | the typical miss between model and reading; big misses count more |
| R² | the share of the spread in the measured temperatures that the model explains (1 = perfect, 0 = no better than one average value; negative = worse than that average) |
| unseen massifs | each massif predicted from a calibration on the other six, as for a glacier without boreholes |
| grid choice | f_ref picked from 8 fixed values (0.03–1.0) at f_ins 1.017 |
| o4, o5 | the 2 m-air versions of 23 Oct, with three (gentle, north, south) or five (gentle, N, E, S, W) orientation parts |
| fb1 | the first free-air version with five parts, before f_sur |
| option A | today's physics: free air, f_sur 1.4 K fixed, five parts, FRHS, f_adv 1.0 |
| FRHS | the central thermal variant: no T_off in the thermal inputs, ERA5's lapse rate, Herron & Langway meltwater, snow refreezing S |
| f_ref, f_ins | the share of the meltwater that refreezes and heats; the scale of the snow's insulation |
| T_off, c_prec, DDF | each glacier's calibrated temperature offset, precipitation factor and degree-day melt factor |
| `cal_fa`, `cal_2m` | the new recalibrations of all Alpine glaciers, on the free air and on the 2 m air |
| orientation factor | how much more or less a part of a band melts than flat ground in the same weather |
| melt-logging run | a GloGEM run that also writes the thermal melt of every band and month |
| emulator | a fast statistical stand-in for GloGEM, trained on the design runs |
| rules commit | the commit that fixes every rule before any fit, so results cannot bend the rules |
| zone | 3000–3800 m, the change from firn to no firn |
| profile level, shape | a hole's mean temperature in one year; how its temperature changes with depth around that mean |
| cold content | the heat it takes to warm the winter snow and firn to 0 °C; refreezing meltwater supplies it |
| permeable depth | how deep meltwater can sink into the firn before it refreezes or meets ice |
| floor | about 1.7 °C: the closest a band model can get to single boreholes (chapter 11) |
