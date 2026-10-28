# Plan: the mass balance on the free air, Alps only: what it takes and how long (28 Oct 2026)

## Short answer

- **About 2 weeks of calendar time, and less than one day of computer time.**
- The cluster is not the bottleneck: recalibrating all Alpine glaciers takes under an hour.
- The time goes into the work around it:
  - a decision on glacier tongues below ERA5's ground (cold pools);
  - the projections;
  - the reviews;
  - redoing the ice-temperature calibration on the new mass balance.
- **Without projections (past only, 1940–2025): about 1.5 weeks.**
- **Free air is an option (your comment, 28 Oct).** The default stays today's 2 m air. Every change below runs only
  when a user switches free air on in the config.

| step | what | working time | computer time |
|---|---|---|---|
| 0 | Where does the air change? (maps, no runs) | ½ day | minutes |
| 1 | The switch in GloGEM (an option, off by default) + max-effort review | 1½–2 days | test runs |
| 2 | Glacier tongues below ERA5's ground (your choice) | 1–1½ days | none |
| 3 | The air's own effect, then the recalibration of the 4,097 Alpine glaciers | ½–1 day | ~1–1.5 h |
| 4 | Check against GLAMOS (seasonal balances) | 1 day | none |
| 5 | Projections (only if needed now) | 1½–2 days | 1–1.5 h |
| 6 | Ice-temperature calibration redone on the new mass balance | 1 day | 4–6 h |
| 7 | Book chapter, decision log, model-changes doc, main | ½ day | — |
| N | The notebook "What the free air does to GloGEM" (grows with each step) | 1–1½ days | none |
| | **total** | **≈ 8–10 working days** | **< 1 day** |

Steps 0 and 2 need no code, so they can run while Step 1 is being built. The notebook is filled step by step.

```
 working day   1    2    3    4    5    6    7    8    9    10
 Step 0        ■
 Step 1        ■■■■■■■■■ (review)
 Step 2             ■■■■■■ (your choice)
 Step 3                        ■■
 Step 4                          ■■■■
 Step 5                               ■■■■■■■ (if needed now)
 Step 6                                      ■■■■ (cluster ~5 h)
 Step 7                                           ■■
 notebook      ■ ·  ·  ·  ·  · ■■ · ■■ · · · ■■ · ■■ ■
```

## Context

- **The thermal model already runs on the free air.** The mass balance still uses the 2 m air plus each glacier's
  calibrated temperature offset (T_off). One model, two airs.
- **The standing rule (23 Oct)** says input is chosen by its match to its own observations. The free air won that
  test at the high stations (+0.6 °C in winter against −5.2 °C), so the mass balance should follow.
- **A hint that the 2 m air shaped today's calibration:** T_off grows with height. Above 3400 m the median is +3.1 °C,
  and 93 % of those glaciers sit at both other bounds (memory, 12 Oct).
  - But the 2 m air's error is a winter error: in summer both airs agree at the high stations. Melt happens in
    summer, so how much T_off really moves is for Step 3 to measure, not to predict.
- **Scope: the Alps (CentralEurope) only** (your answer, 28 Oct). The free air has been checked against stations only
  here.

## What changes in the model (read from the code)

- **One line makes the mass balance's band air:** `glogem.pro:531`
  `tg = temp[...] + (elev − hclim)·dtdz[m−1] + t_offset`.
  - In words: the 2 m air at the climate grid point, plus the lapse rate times the height difference, plus T_off.
  - It feeds the rain/snow split (`accumulation.pro:26-30`), melt (`meltmodel.pro:22-58`), refreezing and debris
    melt.
- **Two more places compute mass balance from the air:** `firnice_reference_forcing.pro:89` (the reference balance
  used by the thermal model) and `initialise_firnicetemp_spinup.pro:133` (firn accumulation in the spin-up).
- **The free-air pieces already exist** in the `freeair-thermal` branch: `functions/freeair_at.pro` and
  `procedures/initialise/setup_firnice_freeair.pro`. Two things stop them being reused as they are:
  - the setup runs only when the thermal model is on, and the calibration switches the thermal model off
    (`settings.pro:526`);
  - the free-air table ends in 2025, so projections need their own anchoring (Step 5).
- **No lapse rate for the air any more.** The free air brings its own height profile, month by month, including
  winter inversions. Today's air is one grid-point value carried up with a fixed 12-month lapse rate (`dtdz`). With
  free air, `dtdz` is used for the air only below ERA5's ground, and only if Step 2 chooses so.
- **Not changed:**
  - the precipitation and its height gradient;
  - the spread of daily temperatures within each month (still from the 2 m air file);
  - the flow model, which only sees the resulting mass balance.

## The principle: free air is an option, off by default (your comment, 28 Oct)

- **One switch in the config:** `air_temperature = 'era5_2m'` (default) or `'freeair'`, plus the path of the
  free-air table. The default is set in `settings.pro`; you set the switch in your config yourself.
- **Everything new is behind that switch:**
  - the mass balance's band air;
  - the calibration;
  - the reference balance and spin-up of the thermal model;
  - the treatment below ERA5's ground (Step 2);
  - the projection anchoring (Step 5).
- **With free air on, the thermal model follows automatically.**
  - The existing thermal-only switch (`firnice_thermal_air`) stays. It reproduces today's Step 5 setup (thermal model
    on free air, mass balance on 2 m air) and is also off by default.
- **No silent fallback.** A region without a free-air table stops with a clear message.
- **Proof that the default is untouched:** with the switch off, byte-identical results in four tests:
  - a past run;
  - a calibration;
  - a GCM projection;
  - a thermal run.

## Step 0: where does the air change? (no runs, half a day)

- For every Alpine glacier band: free air minus today's band air without T_off, per month, averaged over 2000–2020.
- Shown by height and month, with the band area. This answers where and when the mass balance's air moves.
- **Count the glacier area below ERA5's ground** (the cold-pool zone of Step 2).
- One figure for the book.

## Step 1: the switch (1½ days, including a max-effort review)

- **The switch `air_temperature`** as described in "The principle" above.
- **The free-air setup moves out of the thermal block**, so it also runs in calibration. It runs only when the switch
  is on. Settings checks are updated (`settings.pro:584-587` today requires the thermal model to be on).
- **The hook in three places:** `glogem.pro` after line 531, `firnice_reference_forcing.pro:89` and
  `initialise_firnicetemp_spinup.pro:133`. Band air = free air at the band's height + T_off. `freeair_at.pro` and
  the month index `fa_tix` are reused as they are.
- **Where it is built:** a new worktree, branched from `freeair-thermal` at 61ada09.
  - `GloGEM_wt_freeair` stays at 61ada09: the Bayesian session's chains check that commit before every start.
- **Tests:**
  - switch off: results byte-identical to today in the four tests of "The principle";
  - switch on: the logged band air equals the table + T_off (to 5e-6);
  - a calibration of 3 glaciers runs through.
- **Review:** one independent reviewer at max effort. The thermal switch's review found the one-month-late read.

## Step 2: glacier tongues below ERA5's ground (1–1½ days, including your choice)

- **The problem (chapter 13):** below ERA5's smooth ground the "free air" is ERA5's own extrapolation.
  - At six valley stations 400–1140 m below the ground, it is 2.4–6.5 °C too warm in winter.
  - Today's air there is −2.8 to +2.2 °C off.
  - For the mass balance this matters: the long tongues (Aletsch, Gorner, Mer de Glace) lie below their grid cell,
    and too-warm winter air there means rain instead of snow, and winter melt.
- **Candidates, fixed before any scoring:**
  - (a) free air everywhere;
  - (b) free air above ERA5's ground, today's air below it;
  - (c) a gradual change from (a) to (b) over a fixed depth below the ground.
- **The choice follows the standing rule:** the match to station observations, all months, at valley, middle and
  high stations, never by the model's mass balance. The station series of chapter 13 are reused.
- **You decide** on the scored table.

## Step 3: the air's own effect, then the recalibration (half a day to a day; under an hour of computing)

- **First, the air alone:** one past run (1940–2025) with free air and today's calibrated parameters (T_off, DDF,
  c_prec unchanged). About 10–20 min on 24 cores.
  - It shows what the air alone does to the mass balance, before the calibration absorbs it.
  - It separates the air's effect from the calibration's (notebook part 2).
- **A new results folder.** The calibration file name carries no forcing tag, so writing into today's folder would
  overwrite it.
- **Same method, targets and limits as today:**
  - targets: Hugonnet 2000–2020, glacier by glacier;
  - order: c_prec, then DDF, then T_off;
  - the same regional bounds and the standard tolerance (~0.1 m w.e.).
- **Run:** 24 batches on vierzack03/07, about 3–10 min each. Measured basis: Matthias' Alpine calibration took about
  2 h on one core; ours ran at 0.7 s per glacier.
- **Report:** how T_off, c_prec and DDF move with height; how many glaciers sit at their bounds; whether the
  "T_off grows with height" pattern shrinks.

## Step 4: check against independent mass-balance data (1 day)

- **GLAMOS winter and summer balances**, and the balance with height on the measured glaciers. The geodetic total
  matches by construction, so only the seasons can show a difference.
- **Winter snow at 3000–3500 m**, which is 1.3–2 × GLAMOS today.
- Reported, not a gate (standing rule). It shows where the model needs fixing next, e.g. c_prec or the rain/snow
  threshold.

## Step 5: projections (1½–2 days; only if needed now)

- **Today:** the GMIP4 climate files are already corrected to ERA5's 2 m air at the grid point. Lander or the
  protocol made them; the reference period looks like 1980–2014, to be confirmed with Lander.
- **Method:** future band air = the free air of the reference period at the band's height (per calendar month) +
  the GCM's change at the grid point since that same period.
  - So the warming is the GCM's, and the level and height profile are the free air's.
- **Check:** no jump in band air at the 2025/2026 seam (the GMIP4 seam bug of September was such a jump).
- **Runs:** hindcast + 5 scenarios × 8 GCMs, about 15–20 core-hours, so 1–1.5 h on 24 cores.
  - Either on 03/07, or on 04/06 once GMIP4 v4 ends (Thu 29 Oct ~16:00).
  - The flow spin-ups are rebuilt automatically in the new folder: they are stored per glacier, not per
    calibration, so a new folder is required.

## Step 6: the ice-temperature calibration on the new mass balance (1 day; 4–6 h computing)

- **The Step 5 rules are committed, so the redo is mechanical:**
  - the yaml points to the new calibration folder;
  - 200 design + 45 benchmark runs (about 3 h on 03 + 07), then the second stage, the folds and the fit;
  - check whether the five-part orientation tables use the mass balance's melt; if so, rebuild them first (~40 min).
- **Done with the Bayesian session**, after its present Step 5 result (29 Oct morning).
- **That result stays as the reference on the 2 m mass balance.** It shows how much f_ref and f_ins move when the
  mass balance changes air: a direct test of your caution that f_ref carries the mass-balance calibration.

## Step 7: records and main (half a day)

- Book chapter 15, decision-log rows (Step 2 choice, adoption), a model-changes doc entry, memory.
- Push to GloGEM main after your OK, with the option off by default.

## The notebook: "What the free air does to GloGEM" (your comment, 28 Oct)

- **Where:** `glogemflow_icetemp/notebooks/07_freeair_in_glogem.ipynb`, a Python notebook next to 01–06.
  - It reads saved model output only; no model runs inside.
  - The figures show inline and are not written to disk.
  - The plotting code sits in `scripts/icetemp_analysis/freeair/freeair_mb_figs.py`, so the book's chapter 15 figure
    script uses the same code.
- **Written for other modellers:** each part starts with one plain paragraph saying what it shows and why it matters.
  Each figure caption names its data. A glossary closes the notebook.
- **Parts and figures, filled as the steps land:**

| part | figure | from |
|---|---|---|
| 1. The air | Free air vs 2 m air vs stations, monthly bias by station height (high and valley stations) | ch. 13 data |
| | The free air's own lapse rate by month and height, against GloGEM's fixed 12-month `dtdz`, with winter inversions marked | free-air table |
| | Free air minus today's band air across all Alpine bands, by height and month, area-weighted | Step 0 |
| | Glacier area below ERA5's ground; the Step 2 candidates at the valley stations | Steps 0, 2 |
| 2. The air alone | Change in mass balance with height, winter and summer, with today's parameters | Step 3 (first run) |
| | Snow share of precipitation and melt-season temperature by height, old vs new | Step 3 (first run) |
| 3. The calibration | T_off, c_prec and DDF against glacier height, old vs new; glaciers at their bounds | Step 3 |
| | Map of the T_off change | Step 3 |
| 4. The result | Balance with height on Aletsch, Gorner, Rhone, Mer de Glace and the region, old vs new | Step 3 |
| | Winter and summer balances against GLAMOS; ELA and AAR | Step 4 |
| | Sensitivity: change in balance for +1 °C, old vs new | Step 3 (two short runs) |
| 5. The future (if Step 5 runs) | Alpine volume to 2100 per scenario, old vs new | Step 5 |
| 6. Ice temperature | f_ref and f_ins, old vs new mass balance; borehole fit | Step 6 |

- The sensitivity figure needs two more short past runs (+1 °C on both setups, ~10–20 min each).

## Cluster and coordination

- **GMIP4 v4** holds vierzack04/06 until about Thu 29 Oct 16:00.
- **The Bayesian session** runs Step 5 on 03/07 until about 29 Oct morning.
- Step 3 (24 short jobs) fits on 03 alongside; both sessions get a message before any launch.

## Found on the way (not part of this plan)

- **The per-node climate reader** (NorthAsia, LowLatitudes, Antarctic) reads the lapse rate one month late.
  - `gradient_variability_monthly.pro:29` takes `dtdz[1:12]`. The file starts in January and its lapse column
    repeats every year, so January gets February's lapse rate, and so on.
  - Verified in the code and in a NorthAsia file.
- **Possibly also wrong (needs a run to confirm):**
  - the sub-monthly spread in the same file compares latitude with longitude (`:104-105`);
  - it uses `mtt`, which on that path is an empty string (`initialise_results_files.pro:17`). It may pick the wrong
    grid cell or month.
- **LowLatitudes is in the running GMIP4 v4 rerun.** Telling the GMIP4 session is your call.

## Verification

- Switch off: byte-identical results (past run, calibration, projection, thermal run). Switch on: logged band air =
  table + T_off.
- The notebook runs top to bottom on a fresh kernel without errors.
- Calibration: at least as many glaciers reach their target as today; parameters by height reported.
- Projections: no step at 2025/2026; the hindcast with the switch on equals the calibration's past.
- Ice temperature: all runs OK; 1135 points below 15 m; the same scorer as the current Step 5.

## Glossary

| term | meaning |
|---|---|
| free air | ERA5's temperature on pressure levels, read at the band's height |
| 2 m air | ERA5's air 2 m above its smooth model ground, carried to the band with a lapse rate |
| lapse rate (`dtdz`) | how fast air cools with height; GloGEM uses one fixed value per calendar month |
| inversion | air getting warmer with height, typical of winter valleys |
| ELA, AAR | equilibrium-line altitude; accumulation-area ratio |
| sensitivity | how much the balance changes for 1 °C warmer air |
| ERA5's ground | the smooth terrain height of the ERA5 grid cell (~1900–2400 m in the Alps) |
| cold pool | cold air collecting in valleys in winter, below ERA5's ground |
| T_off | the per-glacier temperature offset from the mass-balance calibration |
| c_prec | the per-glacier precipitation factor |
| DDF | degree-day factor: melt per degree above the melt threshold |
| Hugonnet 2000–2020 | the satellite-measured mass change of every glacier, the calibration target |
| GLAMOS | the Swiss glacier monitoring network (winter and summer balances) |
| GCM | global climate model (8 used in GMIP4) |
| seam | the switch from ERA5 to the GCM in 2025/2026 |
| per-node reader | GloGEM's reader of one climate file per grid point (used by three regions) |
| Step 5 | the Bayesian calibration of the ice-temperature parameters f_ref and f_ins |
| f_ref, f_ins | refreezing heat share; snow and firn insulation |
| max-effort review | an independent check of the code by a separate reviewer at the highest effort |
