# Plan for 15–21 Oct: from the firn-ledger result to the snowline recalibration

*Plan of 2026-10-15 (morning), for your sign-off. It builds on yesterday's results: book chapters 8
and 9, the handover note `plans/2026-10-14_handover.md`, and the overall plan
`plans/2026-10-14_snow_firn_probabilistic_plan.md`, which stands. Nothing has been run or changed
since last night.*

## Context: where we stand

| question | outcome | status |
|---|---|---|
| the forcing | FRH (no T_off, ERA5's lapse rate) central; H a full sensitivity | decided |
| meltwater depth limit in cold firn ('front') | failed twice (two stable states) | closed |
| meltwater refreezing in the winter snow on ice | fixes the false frozen tongues, but loses Tête Rousse's frozen bed | off; judged again on realistic snow |
| the firn ledger (water no deeper than the real firn) | correct (gate passed, IDL = Python); **not adopted** by its rule: zone +1 of the 2 needed, plus the cooling guard | off; judged again after the recalibration |
| zone target D (3000–3500 m) | still not met | time box extended (your rule) |

**What all of yesterday's diagnoses point to.** On small high glaciers, GloGEM's mass balance holds
too much snow for too long. Three things show it:
- the end-of-summer snow line is 80–420 m too low;
- GLAMOS shows a higher ELA and less winter snow;
- the firn model loses the firn 20–30 years earlier.

The calibration target can also be too positive: about +0.25 m w.e. a year for glaciers under
1 km², and +1 at Corvatsch, where the outline does not match. Two of the thermal model's fixes,
the ledger and snow refreezing on ice, are physically sound but cannot pay off on that mass
balance. So the **snowline recalibration (Step 3a) is now the decisive step**, and everything else
waits for it or prepares it.

```
 too much snow / too positive target on small high glaciers   <- fixed by 3a (mass balance)
        |                                       |
   firn too extensive and too thick       Tete Rousse too warm
        |                                       |
   ledger (3b) cannot bite              snow refreezing breaks the frozen bed
        \______________ both judged again after 3a ______________/
```

## A finding from this morning (read-only check): the cooling guard

The six band-decades that failed the firn ledger's cooling guard are not a start-up artefact. They
are firn-line bands with thin or young firn over ice:

| band (RGI7) | what the ledger holds | water reaches: ledger / 'hl' |
|---|---|---|
| Tortin 2855–2885 m (01039) | firn only from 1977–78, after good years: 3–4 m w.e. over ice | 8–9 m / 19–29 m |
| Morteratsch 2875 m (02216), Ortles 3225 m (02821) | reference balance near zero: about 1 m w.e. | 4–5 m / 19–24 m |

There 'hl' soaks 19–29 m deep through what is ice, so the ledger is the more plausible model. The
guard, "keep what 'hl' keeps temperate", cannot tell mature temperate firn from a thin young layer
at the firn line.

**Proposal:** before the rerun after 3a, refine the guard and fix it in advance. It would protect
only bands with a mature firn column (at least 15 layers in the ledger, i.e. 15 years of firn),
where temperate firn is real. This does not change yesterday's verdict: the ledger also failed the
zone condition.

## The steps

### 1. Today (Thu): two preparations for the recalibration

**1a. Is the snow's height profile wrong? (a GloGEM output, and a comparison with GLAMOS bands)**
- *Why:* a c_prec/DDF scan, which is what Alexandra's method does, shifts snow and melt as a whole.
  If the snow increases too steeply with height, no scan can fix the snow line on small high
  glaciers. Two clues: the precipitation gradient (1.5 % per 100 m above the ERA5 cell) and the
  nearly flat GLAMOS winter balances at Corvatsch.
- *What:* GloGEM already computes the band winter balance (`elev_bwb`) but never writes it (finding
  K24). Write it as `bwelev_<id>.dat` beside `belev_<id>.dat`, only when
  `write_mb_elevationbands = 'y'` (default 'n'). That is a small commit on GloGEMflow, and no
  default output changes.
- *Then:* compare the model's band winter balances with GLAMOS's (Corvatsch, Hohsaas, Aletsch, Sex
  Rouge, Tortin, Alphubel...), as gradient with height and level, same years. A Python diagnostic
  in a background agent.
- *Effort:* the GloGEM part in the main session (small, byte-identity check); the comparison by an
  agent (Opus, high).

**1b. The snow-line observations, ready for calibration (Python)**
- Sommer et al. transient snow lines → one file per calibration glacier (date, snow line, σ,
  quality), for the 12 calibration glaciers of 1 km² or more. The model's snow line already exists
  (`Snowline_mon`).
- Which glaciers get which target: snow lines where they exist; GLAMOS winter and annual balances
  where they exist (small glaciers); otherwise the geodetic target only. This is decision 2 below.
- The test set for the port: Hohsaas, Aletsch, Fiescher, Gorner, Morteratsch (snow lines and GLAMOS),
  plus one of her published HMA test glaciers if her branch has the data, to reproduce her numbers.

**1c. Book and memory:** the guard finding and these steps.

### 2. Friday: Alexandra's branch

- Read her branch: the snowline calibration inside the old one-script `glogem.pro`. Map what it
  does onto today's procedures:
  - `calibrate_single_glacier.pro`, the phases, the TSLA reader, the weighted RMSE and spline
    minimum;
  - daily or monthly: her calibration ran daily, our thermal runs are monthly.
- Write a **short port plan for your sign-off**: the new procedure
  (`procedures/calibration/calibrate_snowline.pro`), its switch, what stays byte-identical, the gate,
  and how her published numbers are reproduced.

### 3. Next week: port, recalibrate, rerun

1. **Port and gate** on a GloGEM worktree:
   - off: calibrations byte-identical;
   - on: her test glacier reproduced (snow-line RMSE and parameters close to hers).
   - Main session; one independent reviewer at max. The calibration is core GloGEM, so
     the model-changes page gets an entry.
2. **Diagnostic recalibration of the calibration glaciers**, with the targets of decision 2. Compare
   with GLAMOS and the snow lines, period-matched: winter balance, ELA, snow line in spring and at
   the end of summer.
3. **Your decision:** do the thermal runs use the recalibrated parameters? Documented as a deviation
   from the regional calibration until the mass-balance side adopts one.
4. **Thermal reruns on the recalibrated mass balance**, each against its own baseline by a rule
   fixed in advance:
   - FRH (the new baseline);
   - the firn ledger, with the refined guard if you agree;
   - snow refreezing on ice.
   - Zone target D (≥ 2 more zone entities right than H) is measured on the final combination.

### 4. Then (unchanged): Step 5

```
Alpine ensemble with the adopted model -> T2 Bayesian calibration -> Z2-Z4 once -> v1.0 freeze
(ultracode panel, if you opt in) -> Chessjengletscher scored once -> hindcast and projections
-> probabilistic maps (P(cold/polythermal/temperate), P(frozen bed), cold-ice range, by year)
```

## Decisions for you

| # | question | recommendation |
|---|---|---|
| 1 | refine the ledger's cooling guard (mature firn only, ≥ 15 layers) before its rerun | yes, fixed before the rerun, with this morning's finding as the reason |
| 2 | targets for the recalibration: snow lines where they exist (≥ 1 km²); for small glaciers with GLAMOS, its balances instead of the geodetic target; partial outlines (Corvatsch, Hohsaas, Felskinn) flagged | yes, as a diagnostic; whether it becomes the standard is the mass-balance side's call |
| 3 | 1a: write the band winter balance (output only, off by default) and test the precipitation gradient | yes |
| 4 | the spin-up on the 2003 geometry (the thickness diagnostic, waiting since 13 Oct) | after the recalibration, unless you want it now |
| 5 | Daniel's 2024 plan of action (for the probabilistic product) | tell me where it is when convenient |
| 6 | the firn model's removal rule (published maps) | yours and Matthias's call; the numbers are in chapter 9 |

## Effort per task (as agreed)

| task | how |
|---|---|
| 1a GloGEM output; port of the snowline calibration | main session; one reviewer at max for the port |
| 1a comparison, 1b data preparation, recalibration diagnostics | background agents (Opus, high); I check their numbers |
| run launching and watching | background scripts |
| plans, book, model-changes page | main session |

## Verification

- 1a: with `write_mb_elevationbands = 'n'`, every output is byte-identical; with 'y', the new file's
  winter balance equals `elev_bwb`, and its glacier-wide sum equals the winter balance output.
- 1b: the snow-line files reproduce yesterday's per-glacier scores (wRMSE, bias) from the same data.
- The port and the reruns get their own gates and pre-registered rules (sections 3.1 and 3.4),
  each scorer committed before its run set.
- Book chapter 9, memory and commits as each step lands.
