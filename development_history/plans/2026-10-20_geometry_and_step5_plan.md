# The 2003-geometry diagnosis, then Step 5 with FRHS (20 Oct, evening)

*Short plan for your review. Terms in the glossary.*

## Context

- **Done today:** S is adopted, so FRHS is the central thermal variant. S now sits alone on local GloGEMflow's
  branch `snow-refreeze-flow`, and its gate passes. With every switch off it writes byte-identical output to
  GloGEMflow; with S on, it is identical to the stage-2 S run (373 of 373 files).
- **Your order:** the 2003-geometry diagnosis is the last physics check. After it, Step 5 runs once, on the
  final physics.
- **Step 5** is the Bayesian calibration of f_ref and f_ins against the boreholes (f_adv fixed at 1.0).
  - It uses one point per borehole, year and depth bin.
  - It adds the borehole and site offsets you chose, with sizes taken from the data.
  - Every Step 5 number so far used the old forcing and only 7 glaciers.
- **Found while preparing it: the Step 5 code cannot yet use one point per hole and year.**
  - It looks up readings and model values by depth. Several holes at one depth would all get the first one's
    numbers (`emulator.py:216`, `calibrator.py:284`, `validation.py:456`).
  - It also cannot leave out 01161, would rewrite the shared catchment file, stops runs after 30 min (they
    take 20–36), and crashes when one design point fails.

## Step 0: S onto GloGEMflow (15 min)

- Move local GloGEMflow forward to `snow-refreeze-flow`, after checking with the AMOC session (it runs IDL
  from the main tree). Local only, not pushed.
- Note on the model-changes page that S is adopted (coordinated with the AMOC session).
- The campaign runs from a fixed worktree, so later edits to the main tree cannot change it.

## Step 1: the 2003-geometry diagnosis

**The question:** the thermal spin-up (1942–71 climate) and the run (1971–2025) both use the 2003 ice
thickness. The ice was thicker then. Deep tongue ice keeps its spin-up temperature for centuries, so the
thin 2003 column could explain the warm tongues and Tête Rousse.

```
                   spin-up (1942–71 climate)         run (1971–2025)
  today            on the 2003 ice, thickness H       same column
  test 1b          on thicker ice, H + dh             top dh removed at the start: ice now at depth z
                                                      sat at depth z + dh during the spin-up
```

- **1a. A table and figure in Python (half a day, no runs).** For each borehole band:
  - the modelled 2003 thickness;
  - the drilled depth, for holes that reach the bed;
  - the thickness in the measurement year;
  - each set against the band's misfit at depth under FRHS.

  This shows whether thin model ice lines up with the warm misfits.
- **1b. A thermal-only test (about 1 day, with gate and review).** It needs a diagnostic switch in IDL, off
  by default and byte-identical when off.
  - The spin-up runs on H + dh. dh is the thinning from 1971 to 2003, from GloGEM's Δh curve scaled to the
    glacier's volume change. dh is never negative, and never Σb/0.9, which would ignore the ice flow.
  - At the run start, the top dh is removed.
  - Run: FRHS on the whole set, so the climate points stay as in stage 2. Read per band at depth against
    FRHS, mainly the Grenz/Gorner tongue and Tête Rousse.
  - This is an upper bound: all the thinning comes off at once, and the colder past surface is left out.
- **Update (20 Oct, evening; Janosch's decisions after 1a):**
  - **1a done.** The model ice is half as thick as drilled (median 0.50 at 36 holes that reach the bed), but
    the shortfall does not line up with the warm sites (rank correlation 0.07).
  - **1b skipped.**
  - **New 1c:** the model's mass balance before 2000 is 0.5–1 m w.e./a too positive (Aletsch 1972–2003:
    +0.2 modelled against −14.4 observed). Test 1c:
    - first find whether too much snow or too little melt causes it (MB-only run of the 8 long GLAMOS
      series);
    - then a thermal-only test in which the thermal model's mass-balance input before 2000 is pulled
      towards the observations, read per band at depth against FRHS.
- **What it decides:** only how much the fixed geometry matters. The proper fix is a thermal model in a
  dynamic hindcast from 1940, a project of its own. It is not part of this plan (decision 2).

## Step 2: make Step 5 work with one point per hole and year (1–1.5 days, then a max-effort review)

**The example used throughout:** one band at Colle Gnifetti, three readings at 16 m depth.

```
  hole A, 1991:  −10.9 °C
  hole A, 2003:  −10.4 °C
  hole B, 2018:   −9.9 °C
```

**2a. Each reading compared with its own model value.**
1. The model is read for every reading at its own date. So there are three model values here: one for 1991,
   one for 2003 and one for 2018.
2. Today the code pairs reading and model value by looking for "the value at 16 m". That worked while each
   depth appeared only once per band.
3. With three readings at 16 m, that lookup always finds the first one. All three readings would be
   replaced by the 1991 reading and compared with the 1991 model value. The two later visits would be lost
   without any warning.
4. The fix: give the readings numbers (1, 2, 3) and pair them by number, not by depth.
5. With the old pooling nothing changes. We check this by computing the likelihood before and after the
   change, at 5 parameter sets.

**2b. How much the fit trusts each reading: the "error model".**

1. **Why an error model is needed at all.**
   - The calibration looks for the f_ref and f_ins that make the differences between model and readings
     small.
   - A band model cannot match one borehole exactly. Neighbouring holes alone differ by 1.5–1.9 °C (the
     floor).
   - So the calibration has to be told how big a normal difference is. Otherwise it would bend the
     parameters to chase differences that no parameter can remove.
2. **The second job of the error model: which readings tend to be wrong together.**
   - If two readings are off for the same reason, they are not two independent pieces of evidence. One
     reason would be that their borehole sits in a colder spot than the band average.
   - Counting them as independent makes the fit too sure of itself. Its uncertainty comes out too small,
     and well-sampled places pull too hard.
3. **What the error model knows today.** Readings of the same band share one model error (3.5 °C²). It is
   strongest for readings at the same depth and fades over about 15 m.
4. **Why that breaks with one point per borehole and year.**
   - The three readings sit at the same depth in the same band. Today's error model therefore says they
     share their model error *completely*: the model may be off, but by exactly the same amount for all
     three.
   - Then the only difference it allows between them is the thermistor error, 0.1 °C.
   - Suppose the model gets the warming from 1991 to 2018 wrong by just 0.5 °C, which is normal. The fit
     would read that as five times the allowed difference. It would then move f_ref and f_ins to repair it,
     which is wrong and over-confident.
5. **The new terms: each one names something a group of readings has in common.**

| shared by | what it stands for | in the example | size |
|---|---|---|---|
| the same band (as today) | the model's own error for that band, fading with depth | all three readings | estimated in Step 3 (today 3.5 °C²) |
| the same borehole (new) | the hole's spot in the band: a hollow, a ridge, more or less snow. All its readings are off by about the same amount | A 1991 and A 2003 | about 1.5–1.9 °C, estimated in Step 3 |
| the same borehole and year (new) | one visit: its thermistor string, its calibration, that year's weather | each visit on its own | estimated in Step 3 |
| the same site (new) | all bands of one site share its glacier and climate cell | all of Colle Gnifetti's bands; Tête Rousse's 8 bands | estimated in Step 3 |
| the reading itself | the thermistor's precision, now taken from its own hole (today: the band's first hole) | each reading | 0.05–0.2 °C, from the data |

6. **What this changes in the fit.**
   - Readings with much in common count together as less than separate evidence. Independent readings
     count fully.
   - In the example, the two visits of hole A count less than one visit of A plus one of B. Hole B is a
     second spot; A again is the same spot.
   - Tête Rousse's 8 bands count as one site with 8 bands, not as 8 sites.
   - The model's warming between two visits is judged against a realistic visit-to-visit scatter, not
     against 0.1 °C.
7. **What the error model does not change:** the model runs, f_ref and f_ins, and how the model is read at
   each reading's date. It only sets how strongly each reading pulls on the fit.

```
  the three readings at 16 m          what they share                     how the fit counts them
  hole A 1991  −10.9  ┐               band error        (all three)       less than 3 independent
  hole A 2003  −10.4  ┤ hole A        + borehole error  (the two A's)     readings, more than 1:
  hole B 2018   −9.9  ┘               + visit error     (each on its own) A's two visits count less
                                                                          than two (how much less
                                                                          follows from the sizes)
```

8. **Two practical points.**
   - **What counts as one site** comes from a table we check by hand. Grenz is split into its tongue and
     its upper part (Colle Gnifetti): they are far apart and wrong in opposite directions (tongue too warm,
     summit too cold). The two names for Taconnaz are joined.
   - **The leave-one-out checks** (each glacier system held out in turn) use the same error model. The code
     stops with an error if a term is set where it cannot apply.
- **2c. Run infrastructure:**
  - leave out 01161 in the pipeline;
  - leave the shared catchment file alone;
  - timeout 2 h;
  - a failed design point is dropped, not fatal;
  - training spread over the hosts with free cores, then merged;
  - outputs on scratch, not /home (60 runs × 53 MB).

## Step 3: the sizes of the error terms, from the data (half a day)

- **Estimated once** from the stage-2 FRHS residuals: all depths, per hole and year, each massif at its
  chosen f_ref. Method: REML with a free mean, so the +0.6 °C mean bias is not mistaken for site scatter.
- **A simulation check** refits 200 synthetic sets on the real layout. That shows which sizes the data can
  tell apart; sh2 against ℓ is expected to be a ridge.
- **Fixed and committed before any campaign run,** with sensitivity arms:
  - A0: the estimate;
  - A1: no site term;
  - A2: the site term at the top of its range;
  - A3: the other end of the ℓ/sh2 ridge;
  - A4: estimated without Grenz;
  - a reference arm: today's 3.5 °C² over 15 m.
- Conclusions count only if they hold in A0–A3.

## Step 4: the campaign (about 1 day; about 3 h of computing)

- **A new config** `config/step5_frhs.yaml`; your own yaml stays untouched. Contents:
  - the FRHS switches;
  - CentralEurope_fair without 01161;
  - one point per hole and year, target 'depth';
  - f_ref 0.1–1.0, f_ins 0.3–2.0, f_adv 1.0;
  - explained variance 0.999.
- **Stage 1:** 40 runs, then the emulator and MCMC. **Stage 2:** 20 runs around that result, then the
  emulator and MCMC again.
- **Adoption (the rule of 12 Oct, plus two checks):**
  - the gate: R-hat ≤ 1.05; at least 5 runs close to the result; misfit ratio 0.5–2.
  - **new:** the site offsets stay within their expected range, and no single site carries more than 25 % of
    the misfit;
  - held out by glacier system: the RMSE beats the prior-mean run and is at most T1's + 0.10 °C;
  - the emulator check.
- **One real GloGEM run at the result,** scored like stage 2 (per band at depth, sites, frozen beds) against
  FRHS at f_ref 0.85.
- Reported beside it: leave one massif out.

## Step 5: book chapter 12 and memory (half a day)

## Decisions for you

| # | question | recommendation |
|---|---|---|
| 1 | Test 1b on the whole set (dh for every glacier, read everywhere), not only the 4 named glaciers | yes |
| 2 | If 1b shows a large effect, Step 5 still runs on the 2003 geometry. The fixed geometry is then documented as a known structural error, and the dynamic hindcast becomes its own project | yes |
| 3 | The adoption rule of Step 4: the rule of 12 Oct plus the two site checks | yes |

## Glossary

| term | meaning |
|---|---|
| S, FRHS | meltwater refreezing in the seasonal snow on ice bands; FRHS = the central variant FRH with S |
| f_ref, f_ins, f_adv | share of the meltwater that refreezes and heats; snow insulation scale; advection scale (fixed 1.0) |
| Step 5 | the Bayesian calibration: an emulator trained on real GloGEM runs, then MCMC for f_ref and f_ins |
| emulator | a fast statistical stand-in for GloGEM, trained on the design runs |
| design run, stage 1/2 | a GloGEM run at one (f_ref, f_ins); stage 1 spans the whole box, stage 2 sits around the stage-1 result |
| dh | the ice thickness lost since the spin-up period, per band |
| entity, band, site | site × 10 m band × outline, the comparison unit; site = one borehole location; borehole = one hole |
| s2, ℓ | the band's model error and how far it stays correlated with depth (today 3.5 °C², 15 m) |
| sh2, sy2, ss2 | the new shared errors: per borehole, per borehole and year, per site |
| REML | a way to estimate error sizes from residuals that allows for a free mean |
| T1 | the older calibration target, one number per entity; Step 5 must not do worse held out |
| R-hat | an MCMC convergence check, ≤ 1.05 |
| misfit ratio | Mahalanobis distance per point, which should lie between 0.5 and 2 |
| held out by glacier system | each fold leaves out one glacier system (e.g. all of Monte Rosa) and predicts it |
| gate | the check that a new switch, when off, leaves the output byte-identical |
