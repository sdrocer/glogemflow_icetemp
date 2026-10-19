# The missing heat: where we are (19 Oct, ~17:30)

*Update of the approved plan of this evening (archived as `plans/2026-10-19_missing_heat_plan.md`). Short.*

## Context

- **The question:** the right firn (observed) makes the small zone glaciers too cold, so which missing
  heat source was the wrong firn standing in for?
- **The candidates:**
  - **S:** snow refreezing on ice, in the thermal model;
  - **B:** refreezing on ice in the mass balance, via GloGEM's full refreezing model;
  - **W:** milder winter forcing, from the station lapse rates.

## Done since the plan was approved

**Review (max effort) of S and the refreezing fixes.**
- **S was merged correctly.** Thin snow (below about 0.13 m w.e.) got no refreezing; fixed.
- **S's bed check:** S never reaches the bed on the gate glaciers.

**GloGEM's full refreezing model was broken.**
- **The flaw (also on `main`):** its cold reservoir was evaluated once a year, and an October melt month
  used it up, so spring refreezing in the snow never happened on the tongues.
- **The repair** (`01dbb65`, local GloGEMflow):
  - the reservoir is re-armed after every cold month;
  - partial use is booked;
  - two state bugs are fixed.
- **The result:** refreezing on ice bands 1.6 → 48 mm w.e./a, on firn bands 26 → 65 mm w.e./a
  (23 glaciers, 2000–2025, each on its own recalibration).

**The comparison unit is now per band at depth (your point).**
- **That is how calibration and scoring already work:** each reading at its depth and date, in its
  band.
- **Only my site-mean table and the first stage-1 rule averaged.** The rule was amended before any
  stage-1 result was read (`99ad9da`):
  - every band of Corvatsch and Hohsaas must reach an RMSE of at most 1 °C;
  - no other GLAMOS band's RMSE may grow by more than 0.5 °C;
  - every reading is written for depth profiles.

**Running:** stage-1 part 1, six runs (S, W, S + W with today's and with the observed firn), finishing
about 17:45.

## Two changes from your comments

**1. Per band at depth is the default from now on.** Every comparison of measured and modelled is made
per borehole band, at the readings' depths and dates, as the calibration does. Nothing is averaged
over bands or depths; summaries list bands. This is saved as a standing rule (memory and the book).

**2. Stage 1's pass condition, loosened to what the data allow** (amended before any stage-1 result is
read; part 1 is finishing now and stays unread until then):
- **The key bands:** a fix passes if, with the observed firn and mask, every band of Corvatsch and
  Hohsaas reaches an RMSE of at most **1.7 °C**.
  - That is the floor: neighbouring holes differ by 1.5–1.9 °C, and the development glaciers stand at
    1.7 °C.
  - With the right firn alone these bands are 1.7–4.3 °C too cold, so the fix must bring them back to
    the noise level.
  - The first version asked for 1 °C, often for a single reading; that is below the floor.
- **Do no harm:** no other band of the 9 GLAMOS glaciers may get worse than today by more than
  0.5 °C (RMSE).
  - This one stays tight, because it compares two model versions on the same readings, so measurement
    noise does not enter.
- **Stage 2** keeps the standard rules (temperature classes, A's M1 margins, site counting); no new
  thresholds.

## Still to do, in order

| # | what | why |
|---|---|---|
| 1 | read stage-1 part 1 per band and at depth; a depth-profile figure (measured against modelled, per band) for Corvatsch, Hohsaas, Alphubel, Sex Rouge, Tortin | the first answer for S and W |
| 2 | rebase `firn-ledger` onto the repaired refreezing (only after part 1 has finished, since the runs use that worktree); gate (defaults byte-identical) | B in the thermal runs needs the repair |
| 3 | B's mass balance against GLAMOS (winter balance, ELA, band profiles; `compare_snow.py`), C0 against B | B changes the mass balance itself; it must not get worse there |
| 4 | stage-1 part 2: the B pair (FRH, FRHS, FRHLO, FRHLOS on C0 and on B) | is B's extra refreezing the missing heat? |
| 5 | stage 2 for the fixes that pass (7 f_ref, site counting) | adoption by the rules |
| 6 | the book: chapter 10's perfect-firn table replaced by the per-band one, plus the depth profiles; model-changes entries for the two refreezing fixes (`2cce026`, `01dbb65`) | the record, and the main merge |

## Your call later

- **The full refreezing model is a core GloGEM change.** It is repaired and local, not pushed.
  Whether the thermal runs adopt it depends on step 3 (GLAMOS) and step 4.

## Glossary

| term | meaning |
|---|---|
| S, W, B | snow refreezing on ice (thermal); station lapse rates (thermal surface); GloGEM's full refreezing model (mass balance) |
| C0 | the parametrised-refreezing recalibration of the 23 glaciers, the fair baseline for B |
| per band at depth | each borehole reading below 15 m compared with the model in its 10 m band, at its depth and date |
| RMSE (per band) | root-mean-square model minus measured over that band's readings |
| stage 1 / 2 | a quick look at f_ref 0.85; then the full scored test at the 7 values |
