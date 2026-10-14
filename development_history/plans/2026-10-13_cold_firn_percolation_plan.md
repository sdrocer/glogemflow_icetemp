# Short plan: meltwater in cold firn, re-specified (2026-10-13)

## Why

- **Today ('hl'):** each month, meltwater warms the column layer by layer down to the close-off
  depth: 14–99 m at the high sites.
- **Observations:** in cold firn, meltwater refreezes within about 4 m (Colle Gnifetti: Mattea et
  al. 2021, z_lim = 4 m; Col du Dôme 2012: 4–5 m, Gilbert et al. 2014). Refrozen ice layers block
  deeper percolation, and the rest runs off.
- **The first switch failed** (`firnice_perm_limit = 'cold'`, GloGEM `4a5e962`). Its window
  restarted every month at the first node colder than the melting point, which is almost always
  the 1 m node. In practice it was a fixed 4 m cap below the surface, and it cooled every
  temperate firn band by 1.3–4.8 °C (gate runs, review panel).
- **Temperate firn behaves differently:** it is permeable, and water reaches whatever depth the
  winter cold reached.

## The rule (GloGEM branch `perc-wetting-front`, `8e7e49e`, switch `firnice_perm_limit = 'front'`)

At the first wet month of each hydrological year the season's limit is set from the profile itself.

```
 first wet month of the season, looking down from the surface
   a temperate layer starts within 15 m?   (15 m = depth of zero annual amplitude, below the winter cold)
     yes -> water passes it and reaches 4 m below its bottom          (temperate firn stays temperate)
     no  -> water reaches 4 m below the surface; the rest runs off    (cold firn: Mattea et al. 2021)
```

- It holds for the whole hydrological year. Each month the water warms nodes from the top down to
  the limit, never below the close-off depth.
- The limit comes from the temperatures, not from a stored front, so the spin-up's drift
  correction (which extrapolates the temperatures) carries it consistently.
- The two depths are fixed constants from the literature, not calibrated: 4 m (Mattea et al. 2021;
  Gilbert et al. 2014) and 15 m (the seasonal-layer depth used throughout).
- Where it acts:
  - At cold dry sites the limit stays 4 m below the surface.
  - In temperate firn the water passes as with 'hl'.
  - Where a temperate layer survives the winter below a cold skin, it can deepen by up to 4 m a
    year, the kind of warming Grenzgletscher showed (+6.8 °C in 1991–2008; Hoelzle et al. 2011).

## Gate, before it can enter a run set

1. Defaults ('hl') byte-identical.
2. Bands that 'hl' keeps temperate at 10–20 m stay within 0.1 °C of 'hl'.
3. Wet depth at Colle Gnifetti (01225, 01230 at 4400 m and above) and Col du Dôme (00773 at 4200 m
   and above) stays at most 5 m.
4. Mass-balance files identical.

## Then

- A small run set: the adopted forcing with 'hl' and with the new rule, at 7 f_ref values.
- Decided by rule B of the pre-registration, with 2 m and 8 m as a sensitivity if it is adopted.
- Code is developed in a separate GloGEM worktree while the fair runs use the main tree, and
  merged after they finish.

## Outcome (2026-10-14): the gate fails twice; 'hl' stays

- **The first form** (GloGEM `2fef881`, scan at the first wet month) acted as a 4 m cap everywhere.
  - The post-winter profile never shows a temperate layer within 15 m: winter cooling reaches
    14–19 m by 0.1–2 °C even in firn that 'hl' keeps temperate.
  - The columns also start cold.
  - Result: wet depth 4 m at every site; 53 837 of the 53 982 band-years that 'hl' keeps temperate
    at 10 m are cooled, by up to 7 °C.
- **The second form** (`9a5bc58`) reads the end-of-season profile and limits by node-slice tops,
  so a layer could grow one node per season. Its gate run is identical to the first.
  - From the cold start, the end-of-September profile never has a node within 0.01 °C of the
    melting point (0 of 107 profile bands in 54 years; closest −0.014 °C).
  - Started from the 'hl' state, the rule keeps temperate firn temperate in every band-year.
- **So the rule has two self-sustaining states.** Which one a band ends up in depends on where
  it starts and on a 0.01 °C tolerance: 28 % of bands lie within 0.05 °C of it.
  - That cannot be defended in a projection model.
  - The other fix the review offered was a temperate start. It would make the rule a one-way
    ratchet toward cold, and Taconnaz, the case it was meant for, would stay temperate.
- **Gates 1, 3 and 4 passed both times:**
  - defaults byte-identical (373/373 files);
  - Colle Gnifetti and Col du Dôme at most 4 m;
  - mass-balance files identical.
- **What a 4 m cap everywhere does** (the failed runs, for the record):
  - M1 2.11 → 2.80 °C.
  - In the zone, 12 of 13 cold entities are found, but **0 of 6 polythermal** and 0 of 1
    temperate. Overall, polythermal goes 9 → 5 and temperate 3 → 0.
  - Taconnaz, Corvatsch and Tête Rousse improve; Jungfraufirn, Lys, Fiescher and Alphubel break.
- **By rule B 'hl' stays.** The branch `perc-wetting-front` is kept unmerged. The failed runs are
  `failed_perc_FRHW_r0.85_i1.017_a1.00_w0_cap4m[_9a5bc58]`.
- Any limit on meltwater depth needs a criterion that does not depend on the thermal state it
  changes, for example one taken from the climate of the band.
