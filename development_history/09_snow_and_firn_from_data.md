# 9. Snow and firn from data: what the mass balance hands the thermal model

> **In short.** Chapter 8 traced the misses at 3000–3500 m to too much winter snow and too much
> firn on small high glaciers. The same excess snow decides whether meltwater refreezing in the
> winter snow on ice keeps Tête Rousse's frozen bed.
>
> Letting the precipitation factor fall below its 1.3 bound brings the winter snow to GLAMOS's
> values, and with less snow Tête Rousse's frozen bed comes back. But the ELA moves only 10–50 m
> instead of the 120–200 m GLAMOS shows. With only one glacier-wide target, less snow is simply
> balanced by less melt. A second observation must pin where the snow survives: snow lines
> (Alexandra von der Esch's snowline calibration) and firn maps (Matthias Huss's firn model,
> Janosch's first paper).
>
> *Status (2026-10-14): Step 1 done; snow lines, GLAMOS and the firn model compared (Step 2);
> Alexandra's code arrives as a GloGEM branch by the end of the week.*

## What we believed

- The mass-balance calibration (Huss & Hock 2015, as in GloGEM): one target per glacier, the
  geodetic balance 2000–2020. It sets the precipitation factor (bounded 1.3–2.3 in CentralEurope),
  then the degree-day factor, then the temperature offset.
- The thermal model takes from it the snow (insulation), the firn mask (where meltwater enters) and
  the accumulation (burial).
- We assumed the firn mask was good enough. Chapter 6 compared the ELA with WGMS for a few
  glaciers, but never for small high ones.

## What we measured

### The precipitation bound (Step 1; diagnostic, 4 zone glaciers)

Two calibrations of Tête Rousse (00777), Hohsaas (01112), Altels (01265) and Corvatsch (02103)
were run with GloGEM's own calibration, unchanged.
- **'ref':** bound 1.3, as the regional calibration.
- **'low':** bound 0.8, from a copy of the regional parameter table. The shared file and the code
  stay untouched.

The 'ref' run reproduces the regional calibration: c_prec identical, winter balance within 1 %,
ELA within 10 m. DDF and T_off trade off slightly differently, because on 4 glaciers the cell
means that seed T_off differ.

| glacier | c_prec | winter balance, m w.e. (GLAMOS) | ELA 2010s / 2020s, m (GLAMOS) |
|---|---|---|---|
| Corvatsch | 1.30 → 0.95 | 1.20 → 0.87 (≈ 0.83, 2013–25) | 3179 → 3192 / 3272 → 3283 (≈ 3380, 2013–25) |
| Hohsaas | 1.30 → 0.80 | 1.66 → 1.03 (≈ 0.97) | 3222 → 3247 / 3333 → 3345 (≈ 3460, 2018–25) |
| Altels | 1.30 → 0.80 | 1.55 → 0.94 | 2986 → 3043 / 3095 → 3143 |
| Tête Rousse | 1.30 → 0.80 | 1.63 → 1.00 | 3190 → 3201 / 3185 → 3202 |

```
 one target (the geodetic balance), two knobs that trade off:
   less snow (c_prec down) -> less melt needed (DDF down) -> where they balance (the ELA) barely moves
 so the winter snow comes right, the snow line does not
```

The thermal model on both calibrations (FRH, f_ref 0.85; below 15 m, 18 entities, 163 points):

| | 'ref' | 'low' |
|---|---|---|
| RMSE / bias, FRH | 0.96 / +0.65 °C | 0.85 / +0.51 °C |
| RMSE, with snow refreezing on ice (FRHS) | 1.08 °C | 0.95 °C |
| Tête Rousse hole 64 (32.5 m, measured −2.34 °C): FRH / FRHS | −0.43 / −0.07 °C | −0.91 / −0.39 °C |
| Corvatsch, 4 entities | temperate (bias +0.8 to +2.7 °C) | identical |
| robust class changes | | none |

- **Tête Rousse improves.** With less snow, the snow-refreezing switch keeps its frozen bed: hole 64
  is at −0.39 °C, below the melting point. So the switch failed in chapter 8 because of the excess
  snow, as suspected.
- **Corvatsch does not move at all.** Its firn line rises by only 15 m, so the borehole bands stay
  in firn and meltwater keeps them at 0 °C.
- **Hohsaas gets slightly too cold.** The bias at 3275 m goes from −0.86 to −1.21 °C.

### Snow lines and GLAMOS against today's calibration (Step 2; agent diagnosis)

A GloGEM mass-balance run with full output (monthly snow line, ELA, winter and annual balance,
1990–2025) was compared with two datasets:
- the transient snow lines of Sommer et al. (12 of our glaciers; corrected snow line, QA ≥ 0.67,
  weighted by 1/σ, May–October 2000–2025);
- the GLAMOS series.

Biases are model minus observation, in metres:

| glacier (RGI7, km²) | wRMSE | bias May–Jun | Jul–Aug | Sep–Oct | end-of-summer snow line − model ELA |
|---|---|---|---|---|---|
| Aletsch (02596, 82) | 333 | −208 | +54 | +315 | +48 |
| Gorner/Grenz (01225, 57) | 311 | −283 | +37 | +170 | −21 |
| Fiescher (02611, 31) | 294 | −277 | +106 | +108 | +65 |
| Morteratsch (02216, 16) | 317 | −199 | +84 | +273 | +45 |
| Bossons (00773, 11) | 421 | −252 | +13 | +48 | −238 |
| Lys (01205, 11) | 377 | −407 | +39 | +53 | −82 |
| Taconnaz (00774, 5) | 402 | −456 | −84 | +4 | −196 |
| Monte Rosa east face (01230, 4) | 452 | −281 | −442 | −153 | −415 |
| Hohsaas (01112, 1.8) | 204 | −185 | −67 | +6 | 0 |
| Gruben (01127, 1.2) | 226 | −138 | −126 | −51 | −122 |
| Ortles (02821, 1.2) | 283 | −313 | +12 | −78 | −119 |

- **Spring snow lasts too long on every glacier.** The May–June snow line is 140–460 m too low. The
  model's snow line is the end-of-month state, so the true bias is larger.
- **By the end of summer the glaciers split.**
  - The big valley glaciers are right to slightly high (−21 to +135 m).
  - The small or steep high glaciers keep their snow line, and so their firn line, 80–420 m too
    low.
- **The precipitation bound changes only the winter balance** (Hohsaas May–June −185 → −184 m). The
  calibration keeps the annual target and gives back less snow as less melt.

GLAMOS, glacier-wide, model minus GLAMOS over the same years:

| glacier | years | winter balance, regional / 0.8 bound | annual balance | ELA |
|---|---|---|---|---|
| Aletsch | 1990–2025 | +0.34 | +0.09 | +29 m |
| Corvatsch | 2013–25 | +0.38 / +0.05 | +0.74 / +0.81 | −150 / −138 m |
| Hohsaas | 2018–25 | +0.54 / −0.05 | −0.08 / +0.09 | −126 / −116 m |
| Tortin | 2015–25 | +0.08 | −0.14 | −6 m |
| Sex Rouge | 2012–25 | −0.19 | +0.30 | (one year with an ELA on the glacier) |
| Alphubel | 2016–25 | +0.96 | −0.58 | +21 m |

- **The calibration target itself may be too positive for small glaciers.**
  - Corvatsch's geodetic target (Hugonnet et al., 2000–2019) is −0.25 m w.e. a year, but GLAMOS
    measures about −1.3 over 2013–19, when the model gives −0.6.
  - At 3200–3350 m the model's band balance on Corvatsch and Hohsaas is 1.1–1.6 m w.e. too
    positive, with either calibration.
  - The outlines differ: GLAMOS Corvatsch is 0.17–0.23 km², RGI7 0.45 km². The band comparison at
    equal elevations is not affected by that.
- **Measured or inferred:**
  - measured (seasonal stakes): Corvatsch, Sex Rouge, Aletsch;
  - mostly inferred (GLAMOS method 0): Hohsaas and Alphubel;
  - Tortin's winter balance is inferred.

### The firn model against the thermal model's firn (Step 2; agent diagnosis)

The firn model (`firn_model/firn_change.pro`) was ported to Python. The port agrees with the IDL
original cell by cell (to 5e-7, IDL's output rounding) on all ten glaciers. It was then run on the
GLAMOS distributed mass balances of 9 calibration glaciers, 1970–2025, and averaged to GloGEM's
10 m bands. Firn thickness is converted at an assumed 600 kg/m³.

| zone bands | decade | firn model: firn share / thickness | thermal model: firn share / wet depth |
|---|---|---|---|
| Corvatsch 3295–3365 m | 1970s | 0.50 / 1.2 m | 1.00 / 12.8 m |
| | 1990s–2000s | 0 (no firn 1982–2013) | 1.00 / 12.8 m |
| | 2010s | 0.16 / 1.5 m | 1.00 / 12.8 m |
| Hohsaas 3275–3345 m | 1980s | 1.00 / 15.7 m | 1.00 / 36.5 m |
| | 2000s | 0.09 / 5.0 m | 0.99 / 36.4 m |
| | 2020s | 0 | 0.44 / 36.4 m |
| Aletsch (Jungfraufirn) 3415–3535 m | 2010s–2020s | 0.99–1.00 / 48–51 m | 1.00 / 68 m |

- **Corvatsch and Hohsaas: the diagnosis holds.**
  - The thermal model keeps permeable firn 20–30 years too long: Corvatsch until 2023–25, where
    the firn model has none from 1982; Hohsaas until 2022–23, where the firn model loses it in
    1998–2003.
  - It also soaks it far deeper than the firn could be: 13–36 m against at most about 3–17 m.
- **Jungfraufirn: the firn extent agrees.** Only the depth is wrong: water reaches 68 m against a
  firn base of about 48 m.
- **The meltwater depth limit (Herron & Langway close-off times perm_frac, 75–90 m, fixed per band)
  is deeper than any firn column the firn model produces.** Step 3b (firn thickness from the mass
  balance) targets exactly this.
- **Glacier-wide on Aletsch, the thermal model has less firn than the firn model** since the
  1990s (7 km² in the 2020s). So the excess firn is a small-glacier problem, not a general one.
- **Caveats:**
  - the GLAMOS outlines (about 2003) differ from RGI7: Hohsaas covers 0.19 of 1.79 km², all above
    3215 m. GLAMOS uses the Swiss Glacier Inventory, where Hohsaas is separate from the
    Triftgletscher below it; RGI7 01112 counts both (Janosch). The geodetic target of 01112 is
    therefore for the combined outline. The same may apply to other small Swiss glaciers;
  - the records of Fiescher, Gorner and Morteratsch start in 2005–2008, from zero firn;
  - the measured-period grids differ strongly from the fixed-date ones in some years (Gorner 2025:
    +0.11 against −3.14 m w.e.);
  - Corvatsch's 2022 grid misses the summer;
  - the firn model's thickness is capped by its 20-year rule.

## What was wrong

- **Calibrating to one glacier-wide target cannot fix where snow survives.** The precipitation
  factor sets how much snow falls; the snow line also depends on how the melt and the snow change
  with height. Lowering c_prec trades snow for melt and leaves the ELA nearly in place.
- **On small glaciers, the target may be too positive.** The geodetic balances of glaciers under
  about 1 km² can disagree strongly with GLAMOS (Corvatsch −0.25 against about −1.3 m w.e. a
  year). A too positive target keeps the firn line low, whatever the snow.
- **Possible cause, not yet tested:** the precipitation gradient (1.5 % per 100 m above the ERA5
  cell, the regional table). On small glaciers high above their cell it adds much snow at the top.
  The GLAMOS band winter balances of Corvatsch are nearly flat with height.
- **Spring snow persists too long everywhere**, by 140–460 m in May–June. That changes the
  insulation and the timing of meltwater for every glacier, not only the small ones.

## What is next

- **The firn model** (`glogemflow_icetemp/firn_model/firn_change.pro`, Beer et al. 2026), driven
  by the GLAMOS distributed mass balances, on our 9 Swiss calibration glaciers. Its firn extent and
  thickness by band and year are compared with the thermal firn mask and the meltwater depth.
- **The snow lines:**
  - today's calibration against the Alpine transient snow lines (Sommer et al., ESSD preprint
    essd-2026-35) and the GLAMOS ELAs;
  - then Alexandra's snowline calibration, ported from her GloGEM branch into a procedure on
    GloGEMflow.
- **Firn thickness for the meltwater depth:** a bookkeeping per band with the firn model's rules
  (Step 3b of the plan).

## Sources

- Plan: `plans/2026-10-14_snow_firn_probabilistic_plan.md`, `plans/2026-10-14_plan_additions.md`
- Diagnostic calibrations:
  `/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/zonecal/{ref,low}/`; the
  regional-table copy is in `zonecal/data_lowc/`
- Runs: `data/bayescal/zcal{ref,low}_{FRH,FRHS}_r0.85_i1.017_a1.00[_nodeep]_w0`
- GLAMOS: `/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier/{corvatsch,hohsaas}_fix.dat`
- Literature: Huss & Hock (2015); von der Esch et al. (in review, snowline calibration of GloGEM);
  Sommer et al. (2026, ESSD preprint); Beer et al. (2026, *The Cryosphere* 20, 5453–5473)
