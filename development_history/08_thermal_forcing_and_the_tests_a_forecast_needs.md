# 8. The thermal forcing, the data we were losing, and the tests a forecast needs

> **In short.** Two runs on all 24 calibration glaciers showed that, below 15 m, the model predicts
> glaciers it was not developed on about as well as the development glaciers. The exception is
> cold high firn: Mont Blanc came out 15.6 °C too warm. The cause is the mass balance's
> temperature offset T_off, which enters every thermal input. **T_off is not the air temperature;
> it closes each glacier's mass budget.** It grows with height: 93 % of the high glaciers sit at
> both calibration bounds and still need a median +3.1 °C. And it says nothing about which high
> firn is wet.
>
> Taking it out of the thermal energy inputs (run F) fixed Mont Blanc and Colle Gnifetti but
> starved sunny high firn, so F failed the check we had fixed in advance. That check was not fair:
> one shared refreezing fraction penalises the forcing that removes melt.
>
> We then asked how far the model is from its aim. Two things came out:
> - On the development glaciers it is already as close to single boreholes as a model with one
>   forcing per band can be: neighbouring boreholes differ by about 1.5–1.9 °C.
> - On other glaciers it does not beat a one-number elevation rule.
>
> Four tests independent of the thermal model followed:
> - **Mountain stations** decided the lapse rate (ERA5's own). They also showed that ERA5 is about
>   5 °C too cold in winter at 3300–3600 m.
> - **Sunshine** failed as a band-level predictor.
> - **Warming over time:** the model warms only half as fast as measured.
> - **The 3000–3500 m zone:** too much winter snow keeps firn too long, the heat of former firn
>   lingers, and meltwater soaks 14–99 m deep where the literature says 4–5 m.
>
> On the way we found that the calibration data had been losing about a third of its points below
> 15 m.
>
> **Main lesson: a parameter that closes a budget must not be read as physics, and a test is only
> fair if each option gets its own parameters.**
>
> *Status (2026-10-13): runs R/F/H/FH/P scored (icetemp `61cb4bf`); week-1 tests done (`b0bc5dc`);
> data fixes (`cf140d0`); GloGEM switches for the lapse rate, cold-firn percolation and a bed log
> (`4a5e962`) and the run tooling (`d4af4d7`). The fair run set is pre-registered
> (`plans/2026-10-13_fair_test_preregistration.md`); a review panel checks it before any run.*

## What we believed

- **T_off corrects the downscaled air temperature**, so the ice should feel it too (chapter 6).
  The GloGEM papers treat it the same way: a systematic error of the temperature forcing (Huss &
  Hock 2015; Zekollari et al. 2019, 2024).
- **The calibration data held every usable Alpine borehole.** The loader's rules were written for
  a 10–30 m target and for one outline per glacier name.
- **Below 15 m the remaining errors were scattered.** The supervisor's question was whether the
  model predicts ice temperatures below the seasonal layer well enough. On the development
  glaciers the median error was 1.7–2.0 °C, so it looked close.

## What we measured

### Two runs on all glaciers, and where they fail

Run at the prior mean and at the trial calibration (icetemp `ebaaf43`), the model predicted the
other glaciers below 15 m about as well as the development glaciers, except in cold high firn:

| site (RGI7) | measured 15–40 m | modelled | T_off |
|---|---|---|---|
| Mont Blanc summit (00786) | −14 to −17 °C | 0 °C | +8.4 °C |
| Colle Gnifetti on 01230 | −13.0 to −13.4 °C | −8.5 to −9.7 °C | +2.67 °C (inherited from its cell) |
| Colle Gnifetti on 01225 (same saddle) | −13.0 to −13.4 °C | −12.1 to −12.9 °C | +0.67 °C |
| Taconnaz 3400–3800 m (00774) | −3 to −9 °C | 0 °C down to 110 m | +0.12 °C |

### Is T_off the air temperature? Two tests that involve no thermal parameter

- **Same site, two outlines.** The Colle Gnifetti boreholes lie metres apart on 01225 and 01230 and
  measure the same firn.
  - Without T_off the two outlines' air temperatures agree to 0.03 °C.
  - With it they are 2 °C apart, and the modelled firn 3.5 °C apart.
- **Dry firn as a thermometer.** Where there is almost no melt, firn at 10–20 m is close to the
  mean annual air temperature. At the sites with a large T_off, firn minus air is:
  - +0.25 °C without the offset;
  - −5.3 °C with it, i.e. firn colder than its own air, which is impossible without melt.

  The air temperature was rebuilt in Python from GloGEM's ERA5 files and reproduces the logged
  surface temperature to 0.005 °C (`forcing_check.py`, icetemp `c5fd0c5`).

### What T_off is

```
 the mass-balance calibration, glacier by glacier
   c_prec  1.3 ── 2.3    ── budget closed? stop
   DDF     1.75 ── 4.5   ── budget closed? stop
   T_off   no bound      ── until the budget closes
```

Across the 4073 CentralEurope glaciers:

| glaciers with their own measured target | at both bounds | median T_off | T_off ≥ 2 °C |
|---|---|---|---|
| mean elevation above 3400 m (181) | 93 % | +3.1 °C | 73 % |
| reaching above 4000 m (58, all targets) | — | — | 42 glaciers |
| mean elevation below 3000 m | rare | about 0 | rare |

- At the top of the Alps the model keeps too much snow even with the least snow and the most melt
  the bounds allow, and T_off removes it as melt.
- Two opposite causes hide behind that:
  - Mont Blanc has a measured balance of −0.06 ± 0.18 m w.e./a and firn as cold as its air: its
    snow leaves by wind and sublimation.
  - Lys has firn near 0 °C: real melt, more than the degree-day model gives.
- For the thermal model, T_off ranks wet against dry high boreholes at chance (0.5). Melt with
  T_off ranks them no better than melt without it (0.86 against 0.87 above 3300 m).

### The forcing switches and their runs (plan of 2026-10-12, morning)

GloGEM switches, all defaults = before (`7e22150`, `8d10930`, `f6c1453`):
- `firnice_forcing_toff = 'n'`: the thermal energy inputs without T_off;
- `firnice_hl1980 = 'y'`: the published Herron & Langway (1980) close-off depth;
- `firnice_lh_water = 'y'`: latent heat per metre of water;
- `firnice_perm_limit = 'profile'`: a diagnostic.

Five runs at the prior mean (scored in icetemp `61cb4bf`):
- R reproduced today in all 308 files;
- the mass-balance files stayed unchanged under every switch;
- F fixed Mont Blanc (+15.6 → −0.8 °C) and the Colle Gnifetti gap (+3.8 → +0.3 °C).

But F starved sunny high firn: Lys +1.8 → −6.8 °C, Breithorn −0.1 → −3.6 °C, Fiescher
+1.8 → −2.5 °C. The deep RMSE of the other glaciers without Mont Blanc and Taconnaz went from
2.01 to 3.57 °C, and the check fixed before the runs (at most +0.05 °C) failed.

Two notes on that failure:
- Without T_off the model *orders* the sites better: rank agreement with the measured deep
  temperatures 0.68 → 0.75, and on the other glaciers 0.40 → 0.57. What was wrong was the
  *level*.
- The check ran both forcings at one shared f_ref. That penalises the forcing that removes melt,
  because f_ref is exactly the knob that sets how much heat the water brings. Run P (water
  stopping at the model's own 830 kg/m³ depth) left Taconnaz unchanged.

### How far the model is from its aim

Below 15 m, prior mean, scored set (23 glaciers, 110 entities):

| | result |
|---|---|
| ice | RMSE 1.3 °C |
| development (Step 5a) glaciers | RMSE 1.7 °C |
| other glaciers | 3.9 °C (2.0 without Mont Blanc and Taconnaz) |
| regime, three classes | 0.50 (chance 0.33); polythermal found 6 of 13, on the other glaciers 1 of 5 |

- **The floor.** Boreholes on the same outline, within 60 m and 5 years of each other, differ by
  2.2–2.7 °C (RMS), up to 8–10 °C at Colle Gnifetti and Dôme du Goûter. A band model cannot come
  closer to a single borehole than about 1.5–1.9 °C, and the development glaciers are at 1.7 °C.
- **A simple rule**, fitted leaving one glacier out: "warm below X m" reaches the same two-class
  regime skill as the model (0.67).
- **Where the error sits:** Taconnaz 41 %, Monte Rosa 27 %, Mont Blanc 12 %. By elevation the
  regime is right only 42 % of the time at 3000–3500 m.
- **What our data can tell apart:** resampling the 20 glaciers gives a 90 % range of 0.53–0.77 for
  the regime skill and 1.7–3.4 °C for the RMSE.

### Week 1: four tests independent of the thermal model (icetemp `b0bc5dc`)

| test | rule fixed before | result |
|---|---|---|
| lapse rate: ERA5 rebuilt as GloGEM does against MeteoSwiss monthly means 1981–2020 (Jungfraujoch, Corvatsch; Weissfluhjoch, Säntis) | smaller mean annual bias above 3000 m wins; a gap under 0.3 °C keeps the floor | ERA5's own rate 1.80 °C, floor 2.36 °C: **ERA5's own rate**. Both are about 5 °C too cold in winter (DJF) and close in summer |
| sunshine: clear-sky radiation from DEMs with horizons, 76 high entities | band sunshine must add ≥ 0.05 ranking skill, p < 0.05 | band +0.013 (p 0.22): **fails**. Sunshine at the borehole adds +0.094 (p 0.005), which a band model cannot use |
| warming over time: 200 repeat pairs at 7 sites | pairs fixed before any model value was read | modelled ÷ measured warming 0.46 (today), 0.34 (F), 0.48 (H), 0.35 (FH). Grenzgletscher 4250 m measured +5 to +6.8 °C, model +0.6 to +2.3 °C |
| the 3000–3500 m zone and the bed | causes listed before | snow on 1 May 1.3–2× GLAMOS and model ELA 125–150 m too low (Corvatsch, Hohsaas), so firn persists too long; heat of former firn kept at Tête Rousse; meltwater 14–99 m deep. Bed: of 37 holes, 29 measured frozen; model finds 19 (today), 22 (F); no temperate bed modelled as frozen |

Literature for the next step (notes in the week-1 folder):
- **Meltwater in cold firn stays shallow:** about 4 m at Colle Gnifetti (Mattea et al. 2021) and
  4–5 m at Col du Dôme in 2012 (Gilbert et al. 2014).
- **Taconnaz** stays cold through ice advected from above (Vincent et al. 2015; Faillettaz et al.
  2015).
- **Tête Rousse's** frozen tongue holds its water pocket (Vincent et al. 2012, 2015).

### The data we were losing

An audit placed each of the 225 Alpine glenglat holes in GloGEM's RGI7 glacier masks and compared
the result with the GLIMS-id lookup the loader used:

| what | effect |
|---|---|
| holes matched by GLIMS id only | Chessjengletscher (6 holes, 2024–25) never ran: its holes lie in 01161, whose GLIMS id is not glenglat's |
| an entity took its first hole's outline | three Colle Gnifetti holes in 01230 ran on 01225; the Seserjoch holes (01195) ran on 01225 |
| entities needed a reading at 10–30 m | the deepest profiles were dropped: six Grenzgletscher tongue holes (to 200–303 m) and Gornergletscher (to 359 m) |
| holes without derived covariates next to holes with them | their readings were dropped, although the configuration includes estimated and unassessed profiles: CG08-1/2, CG13-2, CG21-1, CG99-2/4, GG99-12, GG08-5, the 2017 Col du Dôme hole, Tête Rousse 2016–23 |

With the three fixes: 116 → 128 entities, 1103 → 1421 points, below 15 m 651 → 887, and 22 → 24
outlines.

## What was wrong

- **T_off was read as physics.** It is the third calibration factor, unbounded, and it closes the
  budget after the other two hit their bounds. It is large exactly where cold firn lives. Feeding
  it into the firn's energy put melt on dry summits and hid the melt the degree-day model misses
  on sunny slopes.
- **The test of F was not fair.** One shared f_ref favoured the forcing the parameter had been
  used with.
- **The loader lost data silently**, in four ways, none of which showed in a summary.
- **The ERA5 lapse-rate floor** was added in June against a warm bias that T_off partly caused.
  Measured mountain air temperatures favour ERA5's own rate.

## What we changed and why

| change | commit |
|---|---|
| holes on the outline that contains them (masks; at a divide the curator's outline); entities split by outline; entry with a reading below 15 m or at 10–30 m; covariates filled per hole; `DataHandler.left_out` lists every hole no run will see, and why (20 now, all with a reason) | icetemp `cf140d0` |
| `firnice_thermal_lapse` ('reanalysis' = ERA5's rate), `firnice_perm_limit = 'cold'` (water reaches at most 4 m into cold firn), bed log `firnice_bed_<id>.dat`, a guard on the switch values | GloGEM `4a5e962`; gate runs byte-identical with the defaults |
| `--fref`/`--fins`/`--fadv` and run tags that cannot reuse a folder with other settings | icetemp `d4af4d7` |
| the four week-1 tests | icetemp `b0bc5dc` |

## What is still open

- **The fair run set.** Two forcings × two percolations at seven values of f_ref, f_ref chosen
  leaving one massif out, decided by fixed margins on two main measures. The measures are RMSE
  below 15 m with each outline weighted equally, and three-class skill. The pre-registration is
  reviewed first.
- **The winter cold bias** of the downscaled air at height. It is documented; a station-based
  correction is tested only after v1.0.
- **Sunshine** is a stated limitation of a band model.
- **Warming over time** is the stop sign. It is measured again after calibration.
- **Chessjengletscher (01161)** is sealed: run, but scored only at the freeze.
- **The firn line at 3000–3500 m** may need the snowline calibration, with Alexandra and Lander.
- **For the mass-balance side**, through Janosch: the T_off–height finding, and the RGI7
  band-aspect coding. GloGEM's RGI7 band files count 0 = N … 4 = S while the code expects
  1 = N … 5 = S. It is not an RGI problem, and it affects only the unused melt model '3'.

## Sources

- Plans: `plans/2026-10-12_thermal_forcing_plan.md`, `plans/2026-10-12_progress_and_forcing_plan.md`,
  `plans/2026-10-13_fair_test_preregistration.md`
- Runs: `data/bayescal/ce25_{pred,R,F,H,FH,P}_w0`; gates `gate_*`, `gateB_*`, `gate2{R,H,FRH,HC}_w0`
- Scripts: `scripts/icetemp_analysis/{full_set_runs,full_set_score,forcing_check,forcing_runs_score,borehole_outlines,lapse_station_check,exposure_check,time_test,zone3000_bed_check,forcing_fair_test}.py`
- Outputs (scratch): `/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/{thermal_forcing,forcing_check,week1}/`
- Station data: MeteoSwiss open data (data.geo.admin.ch, downloaded 2026-10-13)
- Literature:
  - mass-balance calibration: Huss & Hock (2015); Zekollari et al. (2019, 2024); Rounce et al.
    (2020); Schuster et al. (2023);
  - cold firn and exposure: Suter et al. (2001); Suter & Hoelzle (2002); Gilbert et al. (2014);
    Mattea et al. (2021); Hoelzle et al. (2011);
  - hanging glaciers and hazards: Vincent et al. (2015, 2020); Faillettaz et al. (2015);
  - Tête Rousse: Gilbert et al. (2012).
