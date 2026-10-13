# Pre-registration: the fair run set of week 2 (the forcing)

*Written 2026-10-13, before any week-2 run, and revised the same day after a five-reviewer panel
(the revisions are listed at the end). Follows the plan of 2026-10-12
(`2026-10-12_progress_and_forcing_plan.md`) and the week-1 results (icetemp `b0bc5dc`). It is
committed together with the scorer `scripts/icetemp_analysis/forcing_fair_test.py` before the
first run; anything changed after that is reported as a deviation, with its reason.*

## What is decided with these runs

**The forcing (Z1):** whether the thermal model's energy inputs keep the mass-balance temperature
offset T_off ("H", today) or not ("FRH": without T_off, with ERA5's own lapse rate, as the station
test chose).

Not decided here:
- **The percolation:** the 'cold' switch is re-specified this week and tested in a second, small
  run set after this decision (rule B below).
- **The station-corrected forcing:** reported only; its adoption is decided later, by its own rule.

Already decided by fixed rules:
- ERA5's own lapse rate for the offset-free forcing (station test, week 1);
- no sunshine physics in v1.0;
- Herron & Langway as published and latent heat per metre of water, in every variant (bug fixes).

Why FRH is the physically preferred option: two tests that involve no thermal parameter show that
T_off is not a property of the local climate. One is the same site on two outlines (Colle
Gnifetti); the other is dry firn as a thermometer. T_off also grows with height (93 % of the high
glaciers sit at both calibration bounds, median +3.1 °C).

A stated limitation: the offset-free air is itself too cold at 3300–3600 m, by about 1.8 °C over
the year and 5 °C in winter (Jungfraujoch, Corvatsch, 1981–2020).

## Data and code

- **Data:** the borehole rules of icetemp `cf140d0`.
- **Sealed:** Chessjengletscher (01161) is not run in this set; it is scored once, at the freeze.
- **Code:** GloGEM `8c38120`: the switches of `4a5e962`, the station lapse rates, deepening switch
  and bed log of `b17fcfe`, and the profile-output fix of `8c38120`, each committed after its gate
  runs. Run tooling: icetemp `d4af4d7` plus its follow-up.
- **No reuse:** the week-1 runs are not reused.

## The runs

| variant | T_off in thermal energy | thermal lapse rate | role |
|---|---|---|---|
| H | yes | floor (today) | decision |
| FRH | no | ERA5's own | decision |
| HR | yes | ERA5's own | reported: separates T_off from the lapse rate |
| FRSH | no | monthly station rates (Jungfraujoch, Corvatsch) | reported: the station-corrected forcing |

- All four variants use today's percolation and the bug fixes.
- f_ref ∈ {0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00}, with f_ins 1.017 and f_adv 1.0: 28 runs.
- H and FRH again at all 7 f_ref values with f_ins 0.785 and with 1.25: 28 runs.
- H and FRH at f_ref 0.55 without borehole deepening, the bed the product will show: 2 runs.
- That is 58 runs on catchment `CentralEurope_fair`, at most 8 at a time on vierzack03, IDL
  licences checked first. H at 0.55 runs first, alone, as a smoke run.

## Scored set (fixed; the scorer stops if the data side differs)

- **What counts:** points below 15 m, without Colle Gnifetti below 40 m (`FIT_EXCLUSIONS`) and
  without 01161. Expected: 105 entities, 800 points, 19 outlines, 7 massifs with points.
- **Measured classes:** 83 cold, 16 polythermal, 6 temperate.
- **Zone:** entities with a band centre of 3005–3495 m (13 / 6 / 1). The range 2950–3550 m is
  reported as a sensitivity.
- **Finite values:** a point counts only if every run gives it a finite value. The number dropped
  this way is reported.

## Classification

- Each entity is classified over its points below 15 m (`thermal_structure.classify_point` /
  `classify_profile`), the same way for model and data.
- Each reading's tolerance is its own hole's, but at least 0.0002 °C per metre of depth.
- A change between variants is **robust** if it holds at tolerances of 0.05 °C, the hole's own
  value, and 0.5 °C.

## Measures

**Decisive (out of sample, below):**
- **M1:** RMSE below 15 m with each outline weighted equally. Only outlines with at least 3 points
  count; in the zone, only those with at least 3 zone points. M1 also chooses f_ref. M1 over all
  outlines is reported.
- **Class counts:** correctly classified entities in each measured class (n of N), overall and in
  the zone.
- **Warming over time:**
  - the week-1 pairs (`time_test.py` pair list), with holes on the outlines of
    `glenglat_borehole_rgi7.csv` and 01161 left out;
  - at each site, the mean modelled change against the mean measured change;
  - measures: the mean absolute site error, and the ratio of the modelled to the measured means;
  - 90 % intervals from resampling sites.

**Reported, not decisive:**
- M2 (mean recall of the three classes), pooled RMSE and bias, two-class skill, the full class
  table.
- Polythermal entities by layer order (cold over temperate, or the reverse), classified per hole
  profile.
- The time-test slope, with Grenzgletscher 4250 m on its own.
- **The frozen-bed table:**
  - holes that reach the bed, have a known depth, and whose deepest reading is at most 5 m above
    the bottom;
  - the model at that depth and in those months;
  - beds colder than −2 °C kept apart from beds near the melting point;
  - counts by hole and by glacier;
  - beside it, the bed log of the two runs without deepening.
- The glacier type along each glacier, from the bed log: temperate, cold, frozen front with a
  temperate upper part, or the reverse.
- The simple rules (elevation, air temperature), fitted in the same folds.
- Step 5a glaciers and the others shown separately.
- The Colle Gnifetti same-site gap, the dry-firn check, rank agreement, and wet/dry ranking at
  10–20 m.

## Out-of-sample selection

- **Massifs:**
  - Mont Blanc (00773, 00774, 00777, 00786)
  - Monte Rosa (01195, 01205, 01225, 01230)
  - Mischabel–Saas (01112, 01127, 01150)
  - Bernese Alps (01265, 01565, 02596, 02611)
  - Valais west and Diablerets (01039, 01309)
  - Titlis (02669)
  - Écrins (00245)
  - Bernina (02103, 02216)
  - Ortles (02821)
  - Ötztal (03116)

  Only the massifs with scored points form folds (expected 7). The scorer stops on an outline
  missing from the list.
- **Procedure:** for each variant and massif m, f_ref is the grid value with the lowest M1 on the
  other massifs, and m is scored at that f_ref.
  - A tie in M1 goes to the f_ref with the lowest M1 on all massifs.
  - Each time-test site is read at the f_ref chosen for its massif.
  - The whole procedure is repeated at f_ins 0.785 and 1.25.
- **Uncertainty:** paired resampling of massifs, 4000 draws. Each drawn copy is labelled
  separately, and f_ref is chosen again within each draw. Reported: the median and 90 % interval of
  each difference, and the share of draws with the same outcome of A.
- **Robustness:** the whole procedure is repeated with each massif removed in turn.

## Decision rules (fixed)

**A. Forcing.** The two packages, FRH against H, judged on the out-of-sample values of the full set.
FRH is adopted unless one of these holds:
1. **M1:** FRH − H is more than 0.3 °C, overall or in the zone.
2. **Class counts** (net, robust changes only). FRH loses:
   - overall: 2 or more polythermal, 2 or more temperate, or 4 or more entities in all;
   - in the zone: 1 or more of the 7 polythermal or temperate entities, or 2 or more in all.
3. **Warming:** FRH's mean absolute site error exceeds H's, and the 90 % interval of the difference
   lies entirely above 0.

Robustness and its default:
- A is **robust** if its outcome stays the same with each massif removed and at both f_ins values.
- **If A is not robust** (decision of 2026-10-13): FRH becomes the central forcing on the physical
  grounds above, and H is carried through the projections as a full sensitivity.
- Either way, the losing package is kept as a projection sensitivity.

**B. Percolation** (a later run set, after A). The re-specified 'cold' switch replaces 'hl' only if
all three hold:
- it gains, net, at least 2 robust zone entities;
- it breaks none of A's conditions;
- it cools no firn band that 'hl' keeps temperate at 10–20 m to more than 0.2 °C below the melting
  point, where that band lies below the cold-firn limit (3400 m on north faces, 4150 m on south
  faces; Suter et al. 2001).

Otherwise 'hl' stays. If 'cold' is adopted, depths of 2 m and 8 m are run as a sensitivity.

**C. Interaction** (after B). Any variant that beats the choice of A and B on conditions 1–3,
without being worse in the zone, is flagged.

**D. Zone target (v1.0).** v1.0 must classify, net, at least 2 more of the 20 zone entities
correctly than H (robust changes, out of sample). If it does not (decision of 2026-10-13), the
time box is extended for the zone work: percolation, winter snow, and the heat kept from the firn
years.

**Flags:**
- **f_ref at a bound:** 1.00 for H or FRH in 4 or more folds counts as evidence that melt is
  missing (sunshine physics, v2). 0.10 means too much melt heat.
- **A changes when a single outline is removed:** named in the paper.
- **Stop sign:** if the 90 % interval of the warming ratio lies outside 0.5–2, or the sign is
  wrong at most sites, one week of diagnosis comes first.

## Revisions after the panel (2026-10-13), before any run

- **'cold' percolation held back.** As implemented it is a fixed 4 m cap below the surface, and it
  cooled every temperate band in the gate runs. It will be re-specified.
- **M1 restricted to outlines with at least 3 points.** Mont Blanc's 2 points made 75 % of H's M1.
- **M2 replaced by robust class counts as decisive.** The old margin of 0.05 was smaller than one
  entity.
- **The zone defined once and made decisive**, consistently.
- **The time test measured by site means.** The slope is pulled toward zero by the noise between
  holes.
- **Robustness required:** A must hold with each massif removed and at both f_ins values.
- **The bootstrap fixed:** each drawn copy labelled separately, f_ref chosen again per draw.
- **The frozen-bed table** compares the same depths.
- **Variants added:** HR and FRSH (reported), and the runs without deepening.
- **01161** is not run in this set.
- **The reason for FRH** now names the station cold bias.
- **A technical deviation, before any campaign run:** the first smoke run stopped at 01225, whose 30
  borehole bands exceeded IDL's 29 file units. The profile output is now buffered and written once
  per glacier (GloGEM `8c38120`), with every file byte-identical in the gate run. The failed run is
  kept as `failed_fair_H_r0.55_i1.017_a1.00_w0_getlun` and not scored.
