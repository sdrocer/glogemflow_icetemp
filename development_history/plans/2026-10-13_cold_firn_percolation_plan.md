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

## The rule

The wetting front advances at most 4 m into cold firn per melt season, and what it has wetted
stays permeable.

```
 surface ── water ──▶ passes the layers above the wetting front (wetted before)
            z_wf    ◀── wetting front: the deepest node the water has brought to the melting point
            z_wf + 4 m  ◀── this season's limit; cold firn below is reached by conduction only
```

- **Wetting front:** a per-band memory `z_wf`, kept from year to year.
- **Each month:** the water passes and warms nodes down to `z_wf + firnice_perm_cold_depth`, but
  never below the close-off depth. Water that is left over runs off and is logged.
- **After each month:** `z_wf` = the deeper of `z_wf` and the deepest node brought to the melting
  point.
- **The season limit** is fixed at the first wet month of each hydrological year (anchor =
  `z_wf` then), so the front moves at most 4 m per season, not per month.
- **Resets:** `z_wf` = 0 when a band turns from firn to ice; it starts again from the surface when
  the band becomes firn again.
- **Spin-up:** `z_wf` starts at 0 and evolves through the spin-up cycles like the temperatures.

What it does:
- At cold dry sites (Colle Gnifetti, Col du Dôme) the front stays within about 5 m.
- At temperate sites the front reaches the whole wet firn, which stays temperate, as with 'hl'.
- At sites in between, a front can advance 4 m per year. This is the kind of warming
  Grenzgletscher showed (+6.8 °C in 1991–2008; Hoelzle et al. 2011).

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
