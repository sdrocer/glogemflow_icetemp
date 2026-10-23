# Plan: the free air as the base — find why the model gets worse, and fix it with parameters first (23 Oct 2026)

## Context

- **The free air is the better input.** At the 7 high stations, ERA5's free air is +0.6 °C in winter, against
  −5.2 °C for today's 2 m air, and it follows the years better (0.98 against 0.87).
- **The ice-temperature model got worse with it** (chapter 13): 2.71 °C on left-out massifs for free air alone,
  against 1.86. The old cold air had been making up for at least three things:
  - the snow surface's own cooling;
  - the colder past that the deep ice remembers;
  - f_ref pressed to its limit.

**Your decisions (23 Oct):**
1. **The free air becomes the base of the ice-temperature model.** Better-constrained input is kept even though the
   score got worse; the worse score means something in the model needs fixing.
2. **Gilbert's surface cooling stays out** of the model and of all tests.
3. **Fix the causes with the parameter space first, before adding new physics** (your question).

**A standing rule from now on (your comment, 23 Oct):** input data are chosen by how well they match observations of
the input itself, never by whether they improve the model's results. If better input makes the model worse, the model
needs fixing. This holds for any future forcing data too. It goes into memory and the decision log as a rule, not
only for this case.

## Your question: can parameters compensate instead of new physics?

**Yes, if two conditions hold:**
- **(a) The parameter acts where the error is.** What matters is the pattern with depth.
- **(b) The boreholes can tell it apart from the other parameters.** Otherwise the calibration cannot set it.

**What the parameters we have can do:**
- f_ref, f_ins and f_adv change the column with nearly the same depth pattern (correlation 0.9–0.98, measured in
  September). Together they act as **one "warmer/colder" knob**.
- That knob can absorb an error that is about even with depth, such as the free air's +0.5 °C at 15–20 m above
  3800 m. It cannot absorb an error that grows with depth.

**What the deep error needs:**
- The deep error does grow with depth: above 3800 m it rises from +0.5 °C at 15–20 m to +3.8 °C at 100–200 m.
- **The model already has a knob with that pattern: the spin-up's air offset** (`firnice_spinup_dT`). It shifts only
  the climate before the run (1942–71); old deep ice feels it most, the top 20 m hardly.
- As a calibrated parameter it would mean "how much colder the past was". The literature gives about 1.5 °C of
  warming since 1900 in the deep Mont Blanc holes (Gilbert & Vincent 2013).
- **The same idea as the flow model's spin-up (your comment).** There, the 1961–90 mass balance and the deformation
  factor A are perturbed until each glacier matches its inventory volume and length. Here, the spin-up climate is
  made colder until the ice temperatures match.
  - **One difference:** the flow spin-up has a target for every glacier (inventory volume and length). Ice
    temperature has targets only at the 20 borehole sites.
  - So the offset has to be **one value for the Alps** (or a simple rule, e.g. by elevation), calibrated at the
    boreholes and carried to all glaciers, like f_ref.
  - Step 2 reports it per site and massif anyway, to see whether one value is enough.

**The risk is compensation again.** Guards:
- every parameter keeps a physically plausible range;
- a parameter that runs to its bound is read as a missing process, not as a result;
- the values are reported as "effective".

```
 depth      error with free air        one "warmer/colder" knob        past-climate offset
            (above 3800 m)             (f_ref, f_ins, f_adv)           (spin-up air colder)
 15–20 m    +0.5 °C   ■                 shifts  ■■■                     barely  ·
 30–50 m    +2.3      ■■■■              shifts  ■■■                     some    ■■
 50–100 m   +3.3      ■■■■■■            shifts  ■■■                     more    ■■■■
 100–200 m  +3.8      ■■■■■■■           shifts  ■■■                     most    ■■■■■■
```
*(the two knob columns show the expected shape; Step 2 measures it)*

**Five parts per band (added after approval, 23 Oct):** Janosch decided in the other window that v1.0 uses five
parts per band (gentle, north, east, south, west). So the base and the knobs run on five parts, with five-part
factors rebuilt for the free air (`week3/orient_parts5_fa`, orientation session, from the stored free-air pixel
factors) and the hole parts of `week3/orient_parts5`. The comparison on today's air is o5. 130 runs instead of 78.

## Step 1 — the new base (40 runs)

- **Settings:** free air without cooling, on the five-part factors rebuilt for the free air, f_ref at all 8 values,
  parts merged per hole. Tag `fb1`, GloGEM `d783f84` (no change needed).
- **Scored with `freeair_rule.py --tag fb1 --variants FA --base-tag o5 --scheme GNESW`** (reported, not a gate) and
  by depth.
- **Expected:** close to fa1's free air (2.71 °C); the rebuilt factors are the only difference.

## Step 2 — which knob moves which depth (90 runs)

- **Settings:** on the new base, at f_ref 0.40, 0.55 and 0.70 (the values the free air chose), tag `fd2`:
  - spin-up air 1 K and 2 K colder;
  - advection f_adv 0.5 and 1.5;
  - insulation f_ins 0.785 and 1.25.
- **Per knob:** its depth pattern (the change of the error in each depth bin, above and below 3800 m, against the
  base) and its effect per site above 3800 m.
- **The reading, fixed before the runs:** a knob matches the deep error if both hold:
  - its depth pattern above 3800 m correlates with the error's own pattern at 0.9 or more (five depth bins);
  - at its stronger setting it removes at least half of the 30–200 m error, while 15–20 m stays within ±0.5 °C.
- **Also measured:** whether the knobs are separable here (correlation of their depth patterns). That decides
  condition (b).

**Steps 1 and 2 together:** 130 runs in three waves, about four hours unattended, on vierzack03/04/07.

## Step 3 — the parameter space (your decision, after Step 2)

- **What you decide:** which knobs go into the calibration, and with which ranges.
- **My expectation, to be checked by Step 2:** f_ref plus the past-climate offset; f_ins only if Step 2 shows it
  separable; f_adv stays fixed (weak, and the same pattern as f_ref).

## Step 4 — the calibration on the free air

- The Step 5 plan of 20 Oct (error model, sites, leave-one-massif-out benchmark) runs on the free air with the chosen
  parameter space.
- If the past-climate offset is in, the calibration's design and emulator must carry it as a fourth parameter. That
  gets a short revision of the Step 5 plan when we get there.

## Later, not in this plan

- The mass balance on the free air (recalibration, GCM correction).
- The geothermal flux as a further deep knob, if needed (there is no switch for it today).
- A snow-surface energy balance (future work, e.g. with IGM).

## What I build and run

- **icetemp, `scripts/icetemp_analysis/freeair/freeair_diag.py`:**
  - takes the variant (free air with or without cooling), the tags and the f_ref grid as options;
  - reports each knob's depth pattern against the base, the pattern correlation, and the separability.
  - Its dry run must reproduce fa1's free-air numbers before any new run.
- **`week3/chain_freeair_base.sh`:**
  - a copy of `chain_freeair_diag.sh` without `--surface-cool`, tags `fb1` / `fd2`, grid 0.40 / 0.55 / 0.70;
  - merge with the hole parts of `orient_parts_v2` (`merge_part_runs.py --suffix` for the spin-up runs);
  - then the scorers.
- **Records:**
  - decision-log rows: the free air becomes the base and parameters come first; the standing input rule;
  - a feedback memory for the standing input rule;
  - a new book chapter 14 for this work (chapter 13 ends with the test);
  - memory.
  - The plan is archived as `development_history/plans/2026-10-23_freeair_base_plan.md`.
- **Coordination:** I tell the orientation session before committing shared files and before using the hosts.

## Verification

- The dry run reproduces fa1's free-air scores from its existing runs (2.71 °C on left-out massifs).
- Run names match the merge and scorer; all 78 runs report OK; every merged run has 1135 points below 15 m.
- Step 2's depth patterns are computed against the base at the same f_ref, and the reading is applied as written.

## Glossary

| term | meaning |
|---|---|
| free air | ERA5's pressure-level temperature read at the band's height |
| f_ref | the share of meltwater that refreezes and heats |
| f_ins | the scaling of the snow and firn conductivity (insulation) |
| f_adv | the scaling of ice advection (how much ice, and its cold, moves down-glacier) |
| past-climate offset | the spin-up's air shifted colder (`firnice_spinup_dT`); stands for a colder past before the run |
| depth pattern | how a parameter changes the modelled temperature at each depth |
| separable | the boreholes can tell two parameters apart (different depth patterns) |
| effective parameter | a parameter whose value also absorbs missing processes, so it is not a pure physical constant |
| fa1, fb1, fd2 | the first free-air test; the new base; this plan's knob runs |
| orientation parts | gentle, north- and south-facing shares of each band, each with its own melt factor |
| error on left-out massifs | each massif predicted with the f_ref chosen on the others |
| Step 5 | the Bayesian calibration of the thermal parameters (plan of 20 Oct) |
