# Free-air mass balance, Step 2: the air below ERA5's ground (fixed 28 Oct 2026, before any scoring)

Part of `2026-10-28_freeair_mass_balance_plan.md`. Janosch chose the blend depth (400 m) on 28 Oct, before any of
the scores below existed.

## The problem

- Below ERA5's smooth ground (the grid cell's height), the "free air" is ERA5's own extrapolation. At six valley
  stations 405–1139 m below their ground it is 2–6 °C too warm in winter (chapter 13); today's air (2 m air plus
  lapse rate) is closer there.
- At the ground on slopes and cols the free air is right (Bocchetta delle Pisse, at its cell's height: 0.6 °C off
  in winter, today's air 5.4 °C).
- 37.6 km² (1.8 %) of the Alpine glacier area lies below ERA5's ground, two-thirds of it within 200 m of it (Step 0).

## The three options

With `d` = depth below ERA5's ground (`hclim − z`, metres), `FA` = free air, `T2` = today's air
(`T_2m + (z − hclim)·dtdz`):

| option | band air |
|---|---|
| (a) free air everywhere | `FA` |
| (b) step at the ground | `FA` above the ground, `T2` below it |
| (c) gradual change over 400 m | `FA` above the ground; below it `(1 − w)·FA + w·T2`, `w = min(d / 400, 1)` |

T_off is added to all three as before.

## The scoring (by the standing rule: the match to observations of the air itself)

- **Stations:** the 6 valley stations (Samedan, Zermatt, Davos, Ulrichen, Chamonix, Sion), the 3 near the ground
  (Grand St-Bernard, Bocchetta delle Pisse, Aiguilles Rouges) and the 7 above 2800 m. Above their ground all three
  options give the free air, so they differ only at the valley stations; the others are reported to show that.
- **Series:** observed monthly means against both airs at the GloGEM climate point of each station (the free air
  interpolated to the station height; the series of chapter 13), over each station's full overlap.
- **Reported per option and station group:** the mean bias in each season and the year, and the root-mean-square
  error of the monthly values.
- **No automatic winner:** Janosch decides on the table. Not used: any model result (mass balance, ice
  temperature).

## Glossary

| term | meaning |
|---|---|
| ERA5's ground | the smooth terrain height of the ERA5 grid cell (`hclim`, ~1900–2400 m in the Alps) |
| free air | ERA5's temperature on pressure levels, read at the band's height |
| today's air | ERA5's 2 m air at the grid cell, carried to the band with the monthly lapse rate `dtdz` |
| cold pool | cold air collecting in valleys in winter, below ERA5's ground |
| root-mean-square error | the typical size of the monthly differences, large ones weighted more |
