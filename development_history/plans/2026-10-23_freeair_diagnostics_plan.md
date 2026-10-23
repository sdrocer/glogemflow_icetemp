# Plan: what moves the deep ice under the free air (23 Oct 2026)

## Why

The free-air test (chapter 13) failed its rule.
- **What it did well:** above 3800 m it put the top 15–20 m right (−0.5 → +0.1 °C).
- **What went wrong:** the deep ice of the high Mont Blanc sites became 2–4 °C too warm (30–200 m: +2.1 to +3.8 °C).
- **Janosch's question:** could the cause be that the rest of the model was tuned to the old air?

In the test only f_ref was re-chosen. The mass balance stayed on the 2 m air, f_ins and f_adv kept their old values, and
the orientation factors had been built on the 2 m air.

**Approved by Janosch on 23 Oct ("yes lets do it").** This round is diagnostic: reported, no adoption rule.

## What runs

All runs use the free air with Gilbert's surface boundary (FAG, 0.95 × air − 1.4), on the orientation factors
rebuilt for the free air (`week3/orient_pixels_fa`, orientation session). Each run is three parts (gentle, north,
south), merged per hole.

```
                          spin-up air    advection f_adv   insulation f_ins     f_ref            runs
 fa2  (the re-test)       as now         1.0               1.017                all 8 values     24
 fd1  SD1 / SD2           1 K / 2 K colder 1.0             1.017                0.55 0.70 0.85   18
 fd1  ADV05 / ADV15       as now         0.5 / 1.5         1.017                0.55 0.70 0.85   18
 fd1  INS_LO / INS_HI     as now         1.0               0.785 / 1.25         0.55 0.70 0.85   18
                                                                                        total    78
```

**What each setting asks:**
- **fa2:** the test re-run with factors that fit the free air. Scored by the original rule against o4 (`freeair_rule.py
  --tag fa2 --variants FAG`), into its own folder; `fa1` stays the record.
- **SD1, SD2 (colder spin-up):** does the deep ice lack a colder past? In chapter 12 a colder spin-up cooled exactly
  the high sites with long memory.
- **ADV05, ADV15 (advection):** was f_adv tuned to the old air? It moved the deep readings at Bossons, Taconnaz, Grenz
  and Fiescher by up to 1.3–2.0 °C.
- **INS_LO, INS_HI (insulation):** does the snow's insulation, tuned with the old air, set how much winter cold
  reaches the firn?

## How it is read (fixed before the runs)

- **Scorer:** `scripts/icetemp_analysis/freeair/freeair_diag.py`. f_ref is chosen leaving one massif out among 0.55,
  0.70 and 0.85.
- **What it reports per setting:**
  - the error by depth, at and above 3800 m and below;
  - the error on left-out massifs (overall, above 3800 m, zone);
  - the error per site above 3800 m.
- **A setting explains the deep warmth** if, above 3800 m, its mean error at 30–200 m is at most half that of fa2
  (same f_ref grid), while 15–20 m stays within ±0.5 °C.
- **If nothing explains it:** the warmth is not a tuning issue. The free air then waits for a longer spin-up (on the
  long summit records) together with a surface energy balance.

## Order and time

1. **Waves:**
   - wave 1: fa2 and the two spin-up settings (42 runs);
   - wave 2: advection and insulation (36 runs).

   On vierzack03/04/07 (05 is busy), GloGEM `d783f84`, about 40 min per wave.
2. **Merge:** 26 parameter sets, the hole parts of `orient_parts_v2`; the spin-up runs with `--suffix _sd-1.0` /
   `_sd-2.0`, a new option of `merge_part_runs.py`.
3. **Score:** fa2 by the rule; then the diagnostics. About two and a half hours in all, unattended
   (`week3/chain_freeair_diag.sh`).
4. **Then:** results per band at depth to Janosch, chapter 13, and the decision log.

## Glossary

| term | meaning |
|---|---|
| free air, FAG | ERA5's pressure-level temperature at the band's height; FAG adds Gilbert's surface boundary 0.95 × air − 1.4 °C |
| f_ref | the share of meltwater that refreezes and heats |
| f_adv | the scaling of ice advection (how much ice and its cold moves down-glacier) |
| f_ins | the scaling of the snow and firn conductivity (insulation) |
| spin-up dT | the spin-up's air shifted by dT; the reference period is 1942–71 |
| o4, fa1, fa2, fd1 | today's air; the first free-air test (old factors); its re-run on rebuilt factors; this diagnostic round |
| orientation parts | gentle, north- and south-facing shares of each band, each with its own melt factor |
| error on left-out massifs | each massif predicted with the f_ref chosen on the others |
