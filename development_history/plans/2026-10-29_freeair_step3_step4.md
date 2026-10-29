# Free-air mass balance, Steps 3 and 4: the runs and what they report (fixed 29 Oct 2026, before any result)

Part of `2026-10-28_freeair_mass_balance_plan.md`. GloGEM branch `freeair-mb` at 5ea744b (free air with the 400 m
change below ERA5's ground, Step 2). Runs in `mbfreeair/step3/` on vierzack03 (`run_step3.sh`).

## Step 3: the runs

All CentralEurope glaciers (RGI7, 4,097), the same settings in every run; only the air and the calibration differ.

| run | air | calibration | what it shows |
|---|---|---|---|
| `cal_2m` | 2 m air | made here | the reference calibration |
| `cal_fa` | free air | made here | the new calibration (also the one Step 6 uses) |
| `past_2m` | 2 m air | `cal_2m` | today's model |
| `past_fa_air` | free air | `cal_2m` | **the air alone**: the new air on the old parameters |
| `past_fa` | free air | `cal_fa` | the new model |
| `sens_2m`, `sens_fa` | as `past_2m`, `past_fa` | T_off + 1 °C | how much the balance changes for 1 °C warmer air |

- **Why our own 2 m calibration**: on three test glaciers our 2 m calibration differs a little from Matthias'
  (00774: T_off 0.00 against 0.12). A like-for-like reference keeps the air's effect apart from setup differences.
- **Why one job per calibration, not batches**: the calibration's second T_off pass averages the first-pass values
  of the glaciers in each climate cell *of the run*; batches would change those averages.
- **Past runs**: hydrological years 1940–2025, fixed inventory geometry (no retreat, no flow), band balances written
  every year (`write_mb_elevationbands`).

**Reported (Step 3):**
- The parameters: T_off, c_prec and DDF against glacier median height, old against new; the number of glaciers at
  each parameter's bounds; the number that reach their target; whether "T_off grows with height" shrinks.
- The air alone: change in winter and annual balance with height (`past_fa_air` minus `past_2m`), 2000–2020.
- The result: balance with height, old against new (`past_fa` against `past_2m`), for the region and for Aletsch
  (02596), Gorner (01225, the Grenz branch's flowline), Rhone (01706) and Mer de Glace (RGI7 ID looked up by
  outline position before the figure is made).
- Sensitivity: the change of the regional balance for +1 °C, old against new.

## Step 4: the check against GLAMOS (reported, not a gate)

- **Glaciers**: the GLAMOS glaciers matched to RGI7 in `data/snowlines/reference/glamos_vs_geodetic.csv` whose
  GLAMOS area is 0.85–1.15 times the RGI7 outline (the rule of the 16 Oct snow check). Findelen and Adler share
  RGI7 01189: Findelen is used. Glaciers outside the area rule are listed, not counted.
- **Years**: hydrological years 2000–2025. Winter balance only in years with a measured winter balance (GLAMOS
  methods 1, 2, 4, 5); the summer balance is annual minus winter in those years.
- **Compared**: `past_2m` and `past_fa` (and `past_fa_air` for information) against GLAMOS:
  - glacier-wide winter, summer and annual balance: mean of model minus GLAMOS over the glaciers, and the typical
    size of the differences (root-mean-square);
  - winter balance by elevation (GLAMOS's own bands), area-weighted over glaciers and years, with the 3000–3500 m
    zone shown separately (today's model has 1.3–2 times GLAMOS's winter snow there);
  - the equilibrium-line altitude.
- **Why only the seasons can tell**: both calibrations match the same geodetic mass change, so the annual balance
  moves little by construction; the free air changes how the balance splits between winter and summer.
- **No adoption rule**: the result says where the model needs fixing next (e.g. c_prec or the rain/snow threshold),
  not whether the free air is kept. That follows the standing rule (input chosen by its match to its own
  observations; Step 2 and chapter 13).

## Glossary

| term | meaning |
|---|---|
| free air | ERA5's temperature on pressure levels, read at the band's height |
| 2 m air | ERA5's air 2 m above its smooth ground, carried to the band with a monthly lapse rate |
| ERA5's ground | the smooth terrain height of the ERA5 grid cell (`hclim`) |
| T_off | each glacier's calibrated temperature offset |
| c_prec | each glacier's precipitation factor |
| DDF | degree-day factor: melt per degree above the melt threshold |
| winter / summer / annual balance | snow gained over winter; mass lost over summer; their sum over the year |
| GLAMOS | the Swiss glacier monitoring network |
| equilibrium-line altitude | the height where a year's gain and loss balance |
| hydrological year | October to September, named after the year it ends in |
