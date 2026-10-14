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

## Outcome (2026-10-14): not adopted by the rule, by one borehole

The run set is FRHS at 7 f_ref plus one run without deepening, from GloGEM `0f4cd85` (`ab565b7` +
the review's rounding fix). It was scored by `snow_refreeze_rule.py` (icetemp `837c4e9`).

| | FRH | FRHS |
|---|---|---|
| M1 / zone M1 (out of sample) | 2.58 / 1.82 °C | 2.50 / 1.82 °C |
| robust class changes | | none (cold 65 → 62 and zone cold 3 → 2 change with the tolerance) |
| warming: mean site error / ratio | 0.56 / 0.80 | 0.52 / 0.85 (difference −0.13 to 0.00) |
| 1. A's conditions | | **pass** |
| 2. frozen beds found / temperate beds called frozen | 12 of 14 / 2 | 11 of 14 / 1: **fail** |
| 3. fronts of Tête Rousse / Altels frozen (2025) | 1.00 / 1.00 | 1.00 / 0.94: **pass** |
| tongue beds below 2900 m, frozen share: Bossons / Taconnaz | 0.84 / 0.79 | 0.03 / 0.26 |
| Aletsch, Grenzgletscher, Fiescher tongues | temperate | temperate |
| zone target D against H | 0 | 0 |

- **The bed that is lost** is Tête Rousse hole 64, at 32.5 m.
  - Measured −2.34 °C. FRH already gives only −0.43 °C, and the switch takes it to −0.07 °C, past
    the frozen threshold.
  - Two bed holes improve: Grenzgletscher 175 is no longer called frozen, and Alphubel South goes
    from −5.6 to −1.8 °C against a measured −2.1 °C.
- **By the rule, the switch stays off.** The branch `snow-refreeze` is kept unmerged. The false
  frozen tongues of Bossons and Taconnaz remain, as a stated limitation.
- **Why Tête Rousse breaks.** The model gives Tête Rousse far too much winter snow (1.86 m w.e. in
  June; precipitation factor at its 1.3 bound), and on its thin bands that snow fills most of the
  ice column. The zone diagnosis traced the 3000–3500 m misses to the same excess snow. With
  realistic snow, the refreezing heat there would be smaller. So the switch should be judged again
  only once the snow on small high glaciers is right.
