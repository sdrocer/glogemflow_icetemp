# 12. Sun and shade, burial, and the air: three tests, and orientation in the model

> **In short.** The review panel of 20 Oct found that the zone 3000–3800 m cannot be fixed by calibration alone,
> and that Mont Blanc and Monte Rosa want opposite amounts of meltwater heat. We tested three suspected causes, each
> against a rule fixed before its runs:
> - **sun and shade (orientation):** a factor F on the thermal model's meltwater, from an hourly energy balance on
>   the real slopes;
> - **burial:** the model buries the high firn under snow much faster than measured;
> - **the air:** ERA5's air temperature at height.
>
> None of them closes the gap between the massifs. Orientation lowers the overall error on left-out massifs
> (2.51 → 2.28 °C), but removes only about a third of the sun-and-shade pattern in the errors inside sites, so it
> fails its rule, twice. **Janosch kept it in the model**, as three parts per band (gentle, north, south), because
> it is physical and the map must show north against south faces. The burial is 2–20 times too fast at the high
> saddles but moves the borehole temperatures by only 0–0.8 °C. The air is about 2 °C too cold at every high
> station, equally at both massifs; it partly cancels a second error, the missing cooling of the snow surface. Both
> stay known limitations of v1.0.
>
> The build with three parts per band, corrected after a max-effort review, lowers the error on left-out massifs to
> 1.86 °C (above 3800 m from 3.04 to 1.87 °C) and removes the sun-and-shade pattern inside sites; all seven massifs
> choose f_ref 0.85. The zone, the frozen beds and the massif split stay unchanged. West faces get about 12 % more
> meltwater heat than east faces, so five parts (adding east and west) are being tested.
>
> *Status (22 Oct, afternoon): the corrected three-part build is scored (below); the five-part runs are running. Plans: `plans/2026-10-20_step5_plan.md`
> (step 1b, revised 21 Oct, and the new step 1c). Decisions: `decision_log.md`.*

## What we believed

- **The zone's wrong bands have four causes** (review panel, 20 Oct): near the melting point (a matter of
  reporting), too much firn (Corvatsch; fixed by the observed firn), too much meltwater heat on a shaded north side
  (Taconnaz 3685–3785 m), and Tête Rousse too warm at depth.
- **The meltwater heat depends on air temperature only.** A shaded north face melts as much as a sunny saddle. In
  week 1 the sunshine at the hole itself had added skill where the band's sunshine had not (chapter 8).
- **Mont Blanc wants little meltwater heat, Monte Rosa a lot.** In stage 2 the best f_ref was 0.17 at Mont Blanc and
  1.0 (the end of the grid) at Monte Rosa. The best f_ref of a glacier follows the temperature correction T_off of
  its mass-balance calibration (rank correlation 0.78).

## What we measured

### Orientation: an investigation, then Stage 1

**How the factor works.** Only the thermal model's meltwater heat changes; the mass balance stays byte-identical.

```
   air temperature ──► today's melt M ──────────────────────────► mass balance (unchanged)
                            │
                            ▼
   sun on the real slope ─► F = melt on this slope ÷ melt on flat ground   (same weather; table by warmth)
                            │
                            ▼
                     thermal melt = F × M ──► × f_ref ──► heat in the firn and ice
```

- F comes from an hourly energy balance (Giesen & Oerlemans 2012): sun on the slope with its terrain horizon, plus
  heat from the air, over 1981–2020 Jungfraujoch weather.
- On flat ground the energy balance melts exactly what the model's degree-day melt gives in that band (June to
  September). So F carries only slope and shade, not a new melt level.
- F is tabled against the month's warmth (−4 to +6 °C): in a warm summer a shaded slope catches up.
- The investigation (four investigators, a design, two hostile critics, a revision; 20–21 Oct) built the factor and
  predicted from the existing runs that it would remove the shade error inside sites.

**Stage 1** (33 runs, tag `o1`, GloGEM `orientation-factor` `1192935`): F at the borehole bands only, in three
strengths (central, weaker, stronger), at f_ref 0.03–1.00. **Rule B, fixed before the runs** (icetemp `ab5ded4`):
adopt if the shade error inside sites is at least halved, also without Taconnaz and Lys (where the idea came from),
and nothing gets worse overall, as for S.

*The shade error inside sites is shown for every test in the figure at the end of the build section.*

| measure (out of sample, below 15 m) | FRHS | Stage 1, central |
|---|---|---|
| shade error inside sites (°C per +0.1 sunshine) | −1.37 | −0.98 |
| same without Taconnaz and Lys | −1.15 | −1.04 |
| overall error on left-out massifs | 2.51 °C | 2.28 °C |
| without Taconnaz and Lys (also left out of the choice) | 1.36 °C | 1.19 °C |
| above 3800 m | 3.04 °C | 2.73 °C |
| zone 3000–3800 m | 2.51 °C | 2.51 °C |
| regime classes, sites, frozen beds | — | unchanged |
| best f_ref Mont Blanc / Monte Rosa | 0.17 / 1.0 | 0.17 / 1.0 |

- **The rule fails on its first condition.** A reduction of a quarter lies inside the slope's own uncertainty, so it
  cannot be credited to sun and shade.
- **Most of the overall gain is between whole glaciers** (Lys alone gives 85 %), not between the sunny and shaded
  bands of one site. Outside Taconnaz and Lys there is a real but modest gain at Bossons and on the large Monte Rosa
  outline (about 0.2 °C each).
- **The quick estimate from existing runs was a poor guide:** it had predicted the tilt would vanish and the overall
  gain would be smaller.

### Burial: the snow that buries the firn

The thermal model buries its firn at the band's reference-period mass balance (`b_ref`). A throwaway GloGEM branch
printed it at the borehole bands; a literature search (mainly Suter's 2002 thesis, Konrad 2013, Licciulli 2020,
Bohleber 2018, Eichler 2000, Preunkert 2024, Gilbert & Vincent 2013, Vincent 1997/2007, Smiraglia 2000, Mariani
2014, Gabrielli 2016, GLAMOS) gave the measured long-term accumulation at 41 of them.

![Modelled burial against measured accumulation](figures/ch12_burial.svg)

*The model buries most high summits and saddles 2–20 times too fast: Colle Gnifetti (0.14–0.6 measured; 2.6 and,
on its second outline, 6.0 modelled), Dôme du Goûter (0.5 against 2.7), Seserjoch, Lys, Fiescherhorn. It agrees at
the Col du Dôme ice-core site (2.5), the 1994 Grenzgletscher core site (2.7) and Jungfraufirn (2.0).*

**The perfect-burial test** (8 runs, tag `b1`; GloGEM `burial-obs` `abf5cce`; rule icetemp `990cd82`): the measured
accumulation replaced `b_ref` at the 41 bands.
- **The split does not close:** best f_ref Mont Blanc 0.17, Monte Rosa 1.0, as before.
- **Burial moves the borehole temperatures by only 0.0–0.8 °C**, while the massif errors are 1.5–3.7 °C (Mont Blanc
  saddles too warm) and 1.5–7 °C (Monte Rosa saddles and Lys too cold).
- The overall error improves somewhat (2.51 → 2.38 °C; above 3800 m 3.04 → 2.79).

**The re-check of orientation on the measured burial** (25 runs, tag `o2`) fails rule B again: the shade error goes
−1.47 → −1.04. Correct burial makes the tilt slightly larger, so burial was not hiding the sun-and-shade signature.
Together, burial and orientation lower the error on left-out massifs from 2.51 to 2.12 °C.

### The air

The thermal model's air temperature (ERA5 at the glacier's grid point, ERA5's own lapse rate, no T_off) was compared
with measured series at seven mountain stations and three controls (no model runs).

![The model's air minus the measured air, by month](figures/ch12_air_bias.svg)

*The model's air is about 2 °C too cold in the yearly mean at every high station: 4–6 °C in winter, within
±0.6 °C in summer. Monte Rosa (blue) and Mont Blanc (orange) are alike.*

- **Monte Rosa minus Mont Blanc:** 0.0 ± 0.3 °C in the yearly mean, +0.1 ± 0.2 °C in summer. The air does not explain
  the split.
- **It partly cancels a second error.** The model's snow surface is the air temperature capped at 0 °C, with no
  radiative cooling. At Colle Gnifetti the measured air is −12.0 °C (Capanna Margherita, 2003–2018; Mattea et al.
  2021), the measured firn at 15–20 m about −13 °C, the model's air −14.0 °C, and its firn about right. Correcting
  only the air overheats the cold firn, as the station lapse rates (W) did (+2–4 °C at Monte Rosa, chapter 11).
- **This is also why f_ref comes out near 1**, against about a seventh of the melt staying in the firn at Colle
  Gnifetti (Mattea et al. 2021): f_ref partly stands in for the missing surface physics.

### Tête Rousse: the spin-up's memory

The thermal column at Tête Rousse already runs down to the deepest reading (up to 70 m), so a bed at the wrong depth
cannot be the cause of its too-warm deep ice. A spin-up 1 and 2 °C colder than 1942–71 (GloGEM `1d6f970`,
`firnice_spinup_dT`) changes nothing there: every band stays 0.4–2.1 °C too warm, to two decimals. The model holds
its deep ice just below melting (about −0.1 °C) whatever f_ref and spin-up; the measured ice is −0.5 to −2.3 °C.
Elsewhere the colder spin-up cools the cold high sites with long memory (Fiescher −0.6 °C, Mont Blanc summit −0.3,
Bossons −0.2) and leaves the overall error unchanged.

## What was wrong

- **The sun-and-shade pattern in the errors is real** (its 90 % range stays clear of zero), but the melt factor
  explains only a third of it. Something else drives most of it; the candidates are the column flipping to
  water-heated near 0 °C, and less meltwater staying in the firn where melt is large.
- **Taconnaz 3685–3785 m and Tête Rousse share a likely cause:** the meltwater scheme switches a column to
  temperate as soon as water reaches the firn-to-ice depth. It is the one thing all the remaining problem bands
  point to.
- **The massif split is still unexplained:** not sun and shade, not elevation, not burial, not the air.

## What we changed and why

- **Orientation kept, by Janosch's decision against rule B** (21 Oct): it is physical, and the map must show north
  against south faces. The decision log says so; the paper claims the overall gain and only a weak within-site
  confirmation.
- **The build: three parts per band** (step 1c; see below).
- **Known limitations of v1.0, stated and not fixed:** the ERA5 air bias with the missing surface cooling (the
  coupled fix belongs to later work with a full surface energy balance), the over-burial of high saddles, the
  massif split, and the meltwater switch near 0 °C.

## The build: three parts per band

```
        one 10 m band of a glacier                 three runs per glacier            each borehole is scored
  ┌─────────────────────────────────────┐        ┌──────────────────────────┐        against the run of the
  │ gentle (slope < 10°)       F_G ≈ 1  │  ───►  │ run G: every band uses F_G │        part it sits in
  │ north-facing (270–90°)     F_N < 1  │  ───►  │ run N: every band uses F_N │
  │ south-facing (90–270°)     F_S > 1  │  ───►  │ run S: every band uses F_S │
  └─────────────────────────────────────┘        └──────────────────────────┘
```

- **Pixels:** EU-DEM 25 m inside GloGEM's glacier mask (431,000 pixels on the 23 calibration glaciers), with the
  terrain horizons of week 1, joined to the model bands by elevation.
- **The melt each band is matched to:** GloGEM now writes every band's thermal melt and air temperature per month
  (`firnice_write_melth`, off by default); one FRHS run gave them.
- **Checks:** the vectorised sun calculation equals the Stage-1 one exactly; fed the holes' own geometry, the code
  reproduces the 187 Stage-1 factor values to 5·10⁻⁷; the logged band climate equals the Stage-1 rebuild.
- **Holes:** each of the 199 calibration holes takes the part of its EU-DEM pixel (74 north, 64 gentle, 61 south);
  142 of 170 agree with the week-1 fine-DEM planes, most disagreements at the 10° line.
- **Runs:** 24 (three parts × 8 f_ref), tag `o3`; the mass-balance files are byte-identical across the three parts.

![Part factors by elevation](figures/ch12_parts.svg)

*The part factors on two large glaciers. North-facing parts get less meltwater heat with height, south-facing parts
more; gentle parts stay near 1. On the Gorner complex 60 % of the band area faces north.*

### What the build gives

**First build and the review.** The first build (tag `o3`) gave 1.96 °C. A max-effort review (22 Oct) then found that
at the high sites many holes got another factor than intended:
- pixels joined the bands by their EU-DEM height, which lies 8–28 m below GloGEM's own DEM there, so they landed one
  to three bands too low and the holes' own parts were often empty (Col du Dôme got 0.50 instead of about 1.0);
- holes took the part of the single 25 m pixel they sit in, not their week-1 plane (75 m), as the plan says;
- empty parts took the band mean, 0.13–0.15 too low for gentle and south parts.

The corrected build (tag `o4`, icetemp `23effe4`): pixel heights shifted onto GloGEM's DEM cell by cell, hole parts
from the week-1 planes, empty parts from the nearest band's same part; factors, the 24 part runs, merge and scores
redone. Scored like Stage 1 (out of sample, below 15 m, each point from the run of its hole's part); reported, not a
gate, since orientation was adopted on 21 Oct.

![The error on left-out massifs in every test](figures/ch12_results.svg)

| measure | FRHS | Stage 1 (factor at the borehole bands) | first build | **corrected build** |
|---|---|---|---|---|
| error on left-out massifs, all bands | 2.51 °C | 2.28 | 1.96 | **1.86** |
| above 3800 m | 3.04 | 2.73 | 2.12 | **1.87** |
| zone 3000–3800 m | 2.51 | 2.51 | 2.51 | 2.52 |
| without Taconnaz and Lys (also left out of the choice) | 1.36 | **1.19** | 1.34 | 1.26 |
| shade error inside sites (°C per +0.1 sunshine) | −1.28 | −0.90 | −0.25 | +0.16 |
| same, without Taconnaz and Lys | −1.03 | −0.92 | −0.23 | +0.26 |
| scatter of a site's bands about the site mean | 1.63 | 1.63 | 1.72 | 1.72 |
| cold / polythermal / temperate bands right | 61 / 9 / 2 | 61 / 9 / 2 | 62 / 9 / 2 | **63** / 9 / 2 |
| in the zone: cold bands right | 4 of 18 | 4 of 18 | 5 of 18 | 5 of 18 |
| frozen beds found (temperate called frozen) | 11 of 14 (1) | 11 of 14 (1) | 11 of 14 (1) | 11 of 14 (1) |
| best f_ref Mont Blanc / Monte Rosa | 0.17 / 1.0 | 0.17 / 1.0 | 0.23 / 1.0 | 0.23 / 1.0 |
| f_ref chosen in the 7 left-out-massif fits | 0.55–1.0 | 0.70–1.0 | 1.0 in six | **0.85 in all seven** |

*Shade errors with all entities; without the three entities Stage 1 left out of the within-site measures: −1.37,
−0.98, −0.30 and +0.10 (figure below).*

- **The sun-and-shade pattern inside sites is gone,** also without Taconnaz and Lys; it now slightly overshoots
  (+0.16). Each hole sees the factor of its own side of the band.
- **The error on left-out massifs falls from 2.51 to 1.86 °C, above 3800 m from 3.04 to 1.87 °C.** Without Taconnaz
  and Lys the build beats FRHS (1.26 against 1.36), though not Stage 1 (1.19), whose factors were the holes' own.
- **All seven left-out-massif fits choose f_ref 0.85.** The first build had pushed f_ref to its upper limit; that was
  a consequence of its defects, not of the parts.
- **Unchanged:** the zone, the frozen beds, Taconnaz 3685–3785 m and the massif split (0.23 against 1.0).
- **The scatter inside sites grows a little** (1.63 → 1.72 °C), and the north part of a band is on average steeper
  (28°) than the north-facing drill sites (17°): a part's factor is an area mean that a single drill site need not
  match.

![The shade error inside sites in every test](figures/ch12_shade_slope.svg)

*The shade error inside sites (dots; lines: 90 % ranges over sites), without the three entities Stage 1 left out.
Stage 1 and the re-check flatten the tilt by about a third, short of the half rule B asked for; the corrected build
removes it.*

### Five parts instead of three

The three-part split cuts at 270°, right through west, the most common direction of the boreholes: 399 of the 1135
scored points lie on west-facing slopes (8 sites), against 119 north, 62 south and 12 east, and 543 on gentle ones.
A hole facing 260° counted as south, one facing 280° as north. A check without model runs (22 Oct,
`parts_scheme_compare.py`) with five parts (gentle; north, east, south, west as 90° sectors):
- **West really differs from east:** in the same band west faces get a median 12 % more meltwater heat (5–95 %: −3 to
  +49 %), as afternoon sun on warmer air should give. West lies between north and south, east between north and west.
- **The holes' part factors come closer to their own** (correlation 0.81 → 0.83, mean gap 0.116 → 0.108): east
  0.19 → 0.14, south 0.12 → 0.07, west 0.12 → 0.10; north worse, 0.12 → 0.16 (the narrower north sector holds the
  steepest faces). Weighted by readings the gain vanishes (0.073 → 0.075): most readings lie on gentle slopes.

*Running (22 Oct, from 13:31): the 40 runs with five parts (tag `o5`), merged and scored like the three-part build.*

## What is still open

- **Step 5 on the three parts:** every design run becomes three runs (a batch of 40 about 4 h instead of 1.3 h); the
  pipeline merges them per hole before the emulator (`merge_part_outputs`).
- **Three or five parts:** decided on the five-part runs (above).
- **The region-wide factors** for the map (all 4079 glaciers; about 3–4 h of computing plus one GloGEM run).
- **The meltwater switch near 0 °C** (Taconnaz, Tête Rousse) and the massif split: after v1.0.

## Sources

- Plans: `plans/2026-10-20_step5_plan.md` (step 1b revised and step 1c, 21 Oct).
- Investigation: `week3/orientation_investigation/` (workflow result, journal, scripts).
- Runs: `o1` (Stage 1), `b1` (burial), `o2` (re-check), `sd` (spin-up), `mth` (band melt), `o3` (the parts; merged
  `o3_*_oCP`); scores in `week3/orient_stage1/`, `week3/burial_test/`, `week3/orient_recheck/`,
  `week3/orient_parts_score/`; the air check in `week3/airtemp_check/`; the burial literature in
  `week3/burial_research/`.
- GloGEM (branch `orientation-factor`, local): `4f84a94`, `1192935` (the switch), `1d6f970` (spin-up shift),
  `2a3d1fd`, `53af4e5` (band melt output); `burial-obs` `abf5cce` (diagnostic); `diag-bref-print` `b1e51e2`
  (throwaway).
- icetemp: `cc0fa0a`, `7bcbb8f`, `ab5ded4` (Stage 1 tables and rule), `990cd82`, `ab09da1` (burial test), `f799b49`,
  `2039034`, `0d546ac`, `4f7937c` (the build).

## Glossary

| term | meaning |
|---|---|
| FRHS | the central thermal variant: no T_off in the thermal inputs, ERA5's lapse rate, Herron & Langway meltwater, and the snow refreezing S |
| f_ref | the share of meltwater that refreezes and heats; an "effective" parameter |
| orientation factor F | how much more or less a slope melts than flat ground in the same weather; multiplies the thermal model's meltwater only |
| part (G, N, S) | the gentle (slope under 10°), north-facing and south-facing share of a 10 m band; each has its own F and its own model run |
| shade error inside sites | within one borehole site, how much warmer than measured the model is per step of less summer sunshine |
| rule B | the adoption rule fixed before Stage 1: shade error at least halved, also without Taconnaz and Lys, and no harm overall |
| burial (b_ref) | the snow added on top of the firn each year (m water equivalent per year), which carries the firn downward |
| T_off | the temperature correction of each glacier's mass-balance calibration; kept out of the thermal inputs |
| EU-DEM | the European 25 m terrain model |
| left-out massif | the error measure: each massif predicted with the f_ref chosen on the other six |
| zone | 3000–3800 m, the firn–no-firn transition |
