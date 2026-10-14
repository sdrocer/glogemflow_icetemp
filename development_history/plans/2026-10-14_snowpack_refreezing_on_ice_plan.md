# Short plan: the heat of meltwater refreezing in the winter snow on glacier ice (2026-10-14)

## Why

- **The problem.** In the bed log, the tongues of Bossons (00773) and Taconnaz (00774) are frozen
  to the bed at 1750–2900 m, 0.2–3 °C below the melting point. This holds with both forcings and
  at every f_ref and f_adv. These tongues are believed temperate. In the hazard product they would
  show as false "frozen fronts".
- **What the diagnostic runs show** (`week2/diag_*`; chapter 8):
  - The cold is not summit ice carried down: f_ref cools the firn above 3500 m by 5–7 °C, but moves
    the tongue beds by at most 0.06 °C.
  - It is set locally, by a cold mean ice surface (−2.4 to −2.7 °C under the winter snow) over
    thin ice, and by ice from the bands near the snow line.
  - Emergence (upward advection in the ablation area) is in the model and active.
- **What the model leaves out.** On ice bands, meltwater and rain run off without releasing any
  heat (`firnice_temperature_model.pro`, ice branch: `fit_water = 0`). In reality, the first
  meltwater of spring percolates into the cold winter snow on the ice and refreezes, in the snow
  and as superimposed ice on the ice surface. It does so until the snow is at 0 °C, and only then
  runs off. That latent heat warms the snow and the ice surface beneath it every spring. It is a
  textbook control on near-surface temperatures in ablation areas (Cuffey & Paterson 2010, ch. 9).

```
 spring, ice band with winter snow             today                 with this change
   surface  ── melt + rain ──┐                 runs off, no heat      refreezes in the cold snow
   snow (cold) ░░░░░░░░░░░░  │                 stays cold until       and at the ice surface until
   ice surface ─────────────  ▼                conduction warms it    they reach 0 °C; the rest runs off
   ice                                         (frozen tongue)        (ice surface warmer each spring)
```

## The rule (switch `firnice_snow_refreeze = 'y'`, default 'n' = today)

- **Where and when:** on ice bands, in months when the band has seasonal snow.
- **How much water:** the month's melt and rain, times f_ref. This is the same fraction as in the
  firn, so no new parameter is added.
- **How deep:** the water warms the nodes that hold snow, top down, as in the firn loop. It reaches
  down to the node whose slice contains the snow–ice interface, so that node's ice gets the heat of
  the superimposed ice. Whatever is left runs off.
- **Bounds:**
  - Water never passes the interface node: the ice below is impermeable.
  - Nothing happens without snow. In summer the bare ice surface is held at the melting point by
    the surface condition anyway.
- **An upper bound, like 'hl' in firn:** each node is warmed fully to the melting point if the
  water suffices.

## Gate, before any run set

1. Default 'n' byte-identical with the fair campaign (FRH at f_ref 0.85, 373 files).
2. Mass-balance-derived files identical.
3. Firn bands identical in their percolation (wet-depth log unchanged).
4. On the tongues, the ice surface warms in spring, and nowhere does the change cool anything.

## Decision rule (fixed before the run set; scorer committed first)

The run set is FRH with the switch at the 7 f_ref values; f_ref is chosen leaving one massif out,
as in rule A. The switch is adopted if all three hold:

1. **No condition of A breaks against FRH:** M1, robust class losses (overall and zone), warming.
2. **The frozen-bed table does not lose:** no measured frozen bed is lost (12 of 14 found today),
   and at most 2 measured-temperate beds are called frozen (2 today).
3. **The two known frozen fronts stay frozen** (bed types, runs without deepening): Tête Rousse
   (00777; frozen tongue holds the water pocket, Vincent et al. 2012) and Altelsgletscher (01265).

Reported beside the rule:
- the frozen share of the Bossons and Taconnaz tongue beds below 2900 m (the target: below 0.5);
- the Aletsch, Grenzgletscher and Fiescher tongues;
- the zone counts and zone target D.

## What could go wrong

- The warming may also remove the measured cold ice at Tête Rousse and Altels. Condition 3 catches
  that.
- The heat is an upper bound. If it overshoots, the snow and ice surface become temperate too
  early in spring. The test is the borehole points in ice (RMSE 1.3 °C today).
- If the tongues stay frozen, the missing heat is elsewhere: water in crevasses, or water content
  in temperate ice. That is stated, not tuned.
