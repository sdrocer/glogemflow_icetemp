# The missing heat: refreezing in the seasonal snow on ice, and a warmer winter

*Plan of 2026-10-19 (evening), short. Terms in the glossary.*

## Context

- **What batch b19 showed:** with the right firn (observed by GLAMOS), the small zone glaciers come out
  2.3–3.4 °C too cold (Corvatsch +1.8 → −3.4 °C, Hohsaas −0.3 → −2.4). The model's too-extensive
  firn has been supplying meltwater heat the model otherwise lacks.
- **What the model does with meltwater in the seasonal snow (your question):**

  | | firn bands | ice bands with winter snow |
  |---|---|---|
  | thermal model (heat) | yes: snow and firn are one column; the meltwater warms it from the top | **no**: it runs off, no heat (switch S exists, off) |
  | mass balance (mass), default scheme | yes (parametrised) | **no**: its cap is multiplied by the firn mask |
  | mass balance, GloGEM's full refreezing model | yes | **yes**: in the snow's cold content and as superimposed ice (off by default; on GloGEMflow it only runs since today's refreezing fix) |

- **Second candidate:** the surface forcing is about 5 °C too cold in winter at 3300–3600 m (ERA5's
  lapse rate, week 1); the station lapse rates (W) are milder.

## The two fixes for ice bands

- **A. Heat (thermal model): switch S** (`firnice_snow_refreeze`, branch `snow-refreeze`, 14 Oct).
  - The meltwater (× f_ref) refreezes in the cold winter snow and at the ice surface before it runs
    off.
  - It was not adopted then because, under our too-snowy mass balance, it cost Tête Rousse's frozen
    bed.
  - Code: carry its two commits onto `firn-ledger`, checking by hand where they meet the ledger's
    water loop.
- **B. Mass (mass balance): GloGEM's full refreezing model** (`refreezing_full = 'y'`), which already
  refreezes in the seasonal snow on ice.
  - It changes the mass balance, so the 23 calibration glaciers are recalibrated with it (a monthly
    re-fit, about 40 minutes).
  - Its own sub-step still counts 30-day months; it gets the same 365-day fix as the thermal model (a
    core change: model-changes entry).
  - First check the mass balance itself against GLAMOS: winter and annual balance profiles, ELA, and
    how much refreezes on ice bands.

## The test

**Stage 1: a quick look** (f_ref 0.85, about 10 runs, under an hour).

| firn | with S | with W | with S + W | with S + B |
|---|---|---|---|---|
| today's ('hl') | FRHS | FRSH | FRSHS | FRHS-B |
| observed + mask | FRHLOS | FRSHLO | FRSHLOS | FRHLOS-B |

- **Rule, fixed now:** a fix goes on to stage 2 if, with the observed firn,
  - Corvatsch and Hohsaas come within ±1 °C;
  - no other GLAMOS site gets worse by more than 0.5 °C.
- **Also read:**
  - the top row: does the fix alone, with today's firn, make the sites too warm? That would confirm
    the compensation.
  - Tête Rousse's frozen bed.

**Stage 2: the scored test** (passing fixes at the 7 f_ref values, site counting):
1. **The fix itself** against FRH (A's conditions, frozen beds, zone sites).
2. **The ledger with the fix** (the 14 Oct rule).
3. **The perfect-firn reading with the fix** (≥ 2 of the 4 GLAMOS sites). If met, investing in the
   firn lines pays off, together with the fix.

**Effort:**
- the merge and the 365-day change in the full refreezing model get a short review at max effort
  before their gate;
- runs as background scripts;
- results in chapter 10.

## Caveats

- **Corvatsch station data enter the station lapse rates.** These are air temperatures, not
  boreholes, but it is noted.
- **The observed firn exists only for 9 Swiss glaciers.** A fix that is adopted must also work with
  the model's own firn.

## Glossary

| term | meaning |
|---|---|
| S (FRHS) | thermal model: on ice bands the meltwater (× f_ref) refreezes in the seasonal snow and at the ice surface before it runs off |
| B | mass balance: GloGEM's full refreezing model (cold content of 8 layers), which also refreezes on ice bands; needs a recalibration |
| W (FRSH) | the thermal surface with station lapse rates (Jungfraujoch, Corvatsch), milder in winter than ERA5's |
| FRHLO | the firn ledger fed the observed (GLAMOS) balance, with the observed firn mask |
| superimposed ice | meltwater that refreezes on the cold ice surface under the snow |
| compensation | one error hiding another: too much firn supplying heat the model otherwise misses |
| stage 1 / 2 | a quick look at one f_ref value; then the full scored test at the 7 values |
