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
