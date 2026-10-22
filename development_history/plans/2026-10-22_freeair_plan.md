# Plan: test ERA5 free air as the ice-temperature model's air (22 Oct 2026)

## Context

ERA5's 2 m air at the grid point (about 1900–2400 m) is 3.6–5.0 °C colder than ERA5's own free air at the same
height in winter. That is a cold layer the model forms over its smooth snowy ground. GloGEM lapses this value up to
the glaciers, so its air is 4.8 °C too cold in winter at the 7 high stations. ERA5's free air read at the station
height is +0.4 °C in winter, −0.4 °C in summer, and follows the year-to-year ups and downs better (0.98 against 0.87).
Report: `reports/Alpine air temperature beyond ERA5.md`.

**Your decisions (22 Oct):**
- Test the free air in the **ice-temperature model first**. The mass balance follows later if it goes into the model.
- Run **two variants**: free air alone, and free air with a fixed snow-surface cooling from the literature.
- **Test now, adopt by a rule written down before the runs.**
- Use the **full 1940–2025 series**, so the spin-up is tested too. You downloaded it yourself.

This reopens the 21 Oct decision-log entry ("no correction, no sensitivity run"); the log will say so.

```
 today (thermal air)                              new option "freeair"
 ERA5 2 m air at the grid point (~2100 m)          ERA5 temperature on pressure levels
   │  in winter it sits in the model's own            │  = the open air at 1000–5000 m,
   │  cold layer over snow                            │    every month 1940–2025
   ▼                                                  ▼
 + lapse rate × (band height − grid height)        read at each band's own height
   ▼                                                  ▼
 band air: −4.8 °C in winter at the stations       band air: +0.4 °C in winter
                                                      │
                     variant FAC only ───────────────►  surface = 0.95 × air − 3.4 °C, never above 0 °C
                                                        (Gilbert et al. 2014, Col du Dôme; nothing tuned)
 mass balance: unchanged, byte-identical in this test
 precipitation: unchanged (only its rain/snow split in the thermal model follows the new air)
```

**Precipitation stays as it is in this test (your question, 22 Oct).**
- **Nothing to switch to.** The fix works because ERA5 has a second air temperature away from the ground. Its
  precipitation exists only at the model ground; the pressure levels hold no precipitation that falls, so there is
  nothing to read at band height.
- **The amount is already pinned.** The mass-balance calibration (c_prec per glacier) matches it to measured mass
  change.
- **The rain/snow split already moves.** In the thermal model it uses the thermal air, so it shifts with the free air.
  This is part of the test.
- **Burial, the precipitation-linked error we know of, is small for ice temperature.** The perfect-burial test moved
  borehole temperatures by only 0–0.8 °C.
- **Later, with the mass-balance stage:** a short investigation of how precipitation is spread with height (ERA5 and
  GloGEM's gradient against station-based grids such as LAPrec and MeteoSwiss, and the DWD reanalyses). The aim is the
  winter snow at 3000–3500 m (1.3–2 × GLAMOS) and the over-burial of high saddles.

---

## Step 0 — The data (your part is done: both CDS requests finished, 22 Oct)

**What you requested, by the web form** (dataset "ERA5 monthly averaged data on pressure levels from 1940 to
present"; checked against the CDS form definition and its cost check):

| section | setting |
|---|---|
| Product type | Monthly averaged reanalysis |
| Variable | Geopotential, Temperature |
| Pressure level | 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000 hPa |
| Year / Month / Time | request 1: 1940–1983; request 2: 1984–2025 / all / 00:00 |
| Geographical area | North 48.5, West −2, South 42, East 17.5 |
| Data / download format | NetCDF4 / unarchived |

- **Why two requests:** all of 1940–2025 at once costs 148,608 against the CDS limit of 120,000. 1940–1983 costs 76,032 and 1984–2025 costs 72,576.
- **Why these levels:** 400–1000 hPa covers about 100 m to 7 km. The 550 hPa level, which Lander's file lacks, lies
  near 5000 m, just above Capanna Margherita. The finer levels below 750 hPa help at the glacier tongues.

**What I do, after you approve:**
1. **Read the two finished jobs through the CDS API**, with `curl` and the key from
   `/home/jabeer/projects/glogemflow_development/.cdsapirc`. The key is read from the file and never printed. Nothing
   needs installing.
2. **Check each job's stored request against the table above before downloading.** If anything differs, I stop and tell
   you.
3. **Download** to the glogem share, a new folder that follows the existing `monthly/<product>/<region>/` pattern:
   `/itet-stor/jabeer/glogem/climatedata/reanalysis/monthly/era5_pressure_levels/CentralEurope/`
   - `era5_pl_monthly_alps_1940_1983.nc`
   - `era5_pl_monthly_alps_1984_2025.nc`
   - `README.txt`: source, the exact request, download date, licence CC-BY 4.0 with the Copernicus attribution
4. **The key file is protected:** it is readable by you only (`-rw-------`, set 22 Oct) and outside every git repository.
5. **Fallback, only if the API fetch fails:** you download the two files in the browser and copy them into that folder
   under those names.

## Step 1 — I check the data before any model run (main session, about 2 h)

1. **Same numbers as before.** Over 1995–2022 the download must equal Lander's file on the levels both have, apart
   from rounding. This checks that the request is right.
2. **The free-air table.**
   - For every ERA5 grid point GloGEM can use for CentralEurope (the same 0.25° grid as `clim_CentralEurope.mdi`), and
     every month 1940–2025, it holds the free-air temperature at fixed heights: 500 to 5000 m, every 100 m.
   - Each value is interpolated linearly in geopotential height between the two levels around it, the same method as
     the station test.
   - Output, a test product so it stays on scratch until adopted:
     `/scratch-fourth/jabeer/tmp_icetemp_trees/week3/freeair/freeair_CentralEurope.nc`
   - Heights below the ERA5 ground use ERA5's own extrapolated levels; this is flagged and reported, not hidden.
3. **Station rescore with the full series.**
   - The 10 stations of chapter 12 over all their months 1971–2025, now including the 550 hPa level.
   - Reported per station, by season and month.
4. **Is the early free air good enough for the spin-up?** (1942–71, when ERA5 had far fewer observations)
   - Compared against long homogenised summit records: Säntis, Jungfraujoch, Grand St-Bernard (MeteoSwiss open data),
     Sonnblick (GeoSphere) and Zugspitze (DWD).
   - **Criterion, fixed now:** at each station, the free-air bias in 1942–71 lies within ±1 °C of its 1995–2022 bias,
     for the year and for winter.
   - If a station fails, I report it before the runs and we decide together whether the spin-up result can be trusted.
5. **Winter cold pools (report only).**
   - Winter free air against valley stations near glaciers that sit in cold pools, for example Samedan, Zermatt and
     Ulrichen.
   - This shows how much too warm free air would be at low glacier tongues. It matters for the mass-balance stage later.

## Step 2 — The switch in GloGEM, ice-temperature model only (main session; one independent reviewer at max effort)

- **Where:**
  - New branch `freeair-thermal` from `orientation-factor` 53af4e5, in its own worktree `GloGEM_wt_freeair`.
  - `GloGEM_wt_orient` stays untouched; its o4 scoring is running.
- **New settings in `procedures/initialise/settings.pro`**, all off by default, in the `firnice_*` block, with checks
  in Zone 2:
  - `firnice_thermal_air = 'era5_2m'` (default) or `'freeair'`;
  - `firnice_freeair_file = ''`;
  - `firnice_surface_cool = 'n'` (default) or `'gilbert'`.
  - `'freeair'` requires `firnice_forcing_toff = 'n'` and `firnice_thermal_lapse = 'reanalysis'`, so no second shift
    is added on top, plus monthly steps and a file that exists.
- **New procedure `procedures/initialise/setup_firnice_freeair.pro`.**
  - It reads, once per glacier, the profile series of the glacier's own grid point (`clim_ilon` / `clim_ilat`, the
    same point the mass balance uses), matched by exact lon/lat.
  - It aligns that series with GloGEM's monthly time axis (`cyear` / `cmon`).
  - A small function then returns the free air at each band's *current* height for a given month, so a thinning
    glacier gets the right air.
- **Where the thermal air is replaced** when `'freeair'` is on: the same five places the 'stations' variant touches.
  The mass-balance `tg` in `glogem.pro:537` is not touched.
  - `processing/firnice_thermal_forcing.pro:20`
    - `ftw_tm` becomes the free air (plus the spin-up shift `ftw_dsp`, unchanged).
    - It then feeds the rain/snow split, the thermal melt, the logged band air and the orientation lookup.
  - `initialise/firnice_reference_forcing.pro`: the reference-period climatology for the thermal melt and rain, and
    `hl_tk` (`:116-117`) for the Herron–Langway depth.
  - `initialise/setup_firnice_orient.pro:18-28`: the 1981–2020 mean the orientation lookup measures warmth against.
  - `initialise/initialise_firnicetemp_spinup.pro:37-44`: the mean annual air for the starting profile and `hl_tk`.
  - The surface itself: `processing/firnice_temperature_model.pro:283` then sees the free air through `tgs_th`.
- **Surface cooling (FAC only):**
  - At `firnice_temperature_model.pro:283`: surface = min(0, 0.95 × air − 3.4).
  - The same in the starting profile (`initialise_firnicetemp_spinup.pro:198`).
  - **Before coding I read Gilbert et al. 2014 to confirm the exact form** (daily vs monthly means, the cap) and note
    it in the commit.
- **Gate:**
  1. Switch off: byte-identical to 53af4e5 on the gate glaciers (`compare_run_files.py`).
  2. Switch on: all mass-balance files byte-identical.
  3. Harness: the logged band air (`firnice_write_melth`, `tairth`) equals the free-air table at the band heights to
     0.01 °C, and the surface obeys the Gilbert form.
- **Review:** one independent reviewer at max effort reads the IDL change and the gate before any run. Every place the
  air enters must use the same air.
- **icetemp:** `full_set_runs.py` / `runner.py` get `--thermal-air freeair --freeair-file … --surface-cool gilbert`,
  written into the run config and the run name (`_fa`, `_fac`).
- **Starting point:** Step 2 can begin as soon as you approve, in parallel with Steps 0 and 1, using Lander's 1995–2022
  file as a stand-in for the harness only.

## Step 3 — The rule, committed before the runs (icetemp plan + scorer)

**Changed 22 Oct afternoon (your decision): three variants.** Gilbert et al. 2014 (Fig. 7) measured the skin as
0.95 × air − 3.4, but their firn model needed +2.0 K on it for penetrating sunlight (Trad). So:
- FA: free air alone (reported);
- FAS: free air + skin, 0.95 × air − 3.4 (reported);
- **FAG: free air + Gilbert's model boundary, 0.95 × air − 1.4 — the candidate the rule judges.**

Scorer `scripts/icetemp_analysis/freeair/freeair_rule.py`; its dry run on o4 reproduces o4's scores exactly.

**Step 1 results (22 Oct):** the spin-up check passes at all 5 summit records. Cold pools: the free air is 2.4–6.5 °C
too warm in winter at valley stations 400–1140 m below ERA5's ground. Only 1 of 105 borehole bands lies below its grid
cell (Grenzgletscher@2015, 10 points), so this test hardly sees it. It matters for the Alps-wide map and the
mass-balance stage.

- **Baseline:** o4, the current model (FRHS, three parts per band), with the same scorer and the same 1135 points below
  15 m, per band at depth and date.
- **Candidate FAG (free air + surface cooling) is adopted if:**
  1. **Overall:** the error on left-out massifs is not worse than o4 by more than 0.2 °C (thermistor precision).
  2. **Height classes:** the same holds above 3800 m and in the zone 3000–3800 m.
  3. **Classes:** regime classes (cold / polythermal / temperate) and measured frozen beds lose at most one *site* net.
     A single band or bed cannot veto, as for S.
- **FA (free air alone) is reported, not a gate.** It shows how much the surface cooling matters.
- **Also reported:**
  - the f_ref chosen, and whether it leaves the upper limit of 1.0;
  - the Mont Blanc / Monte Rosa split;
  - the Colle Gnifetti and Col du Dôme bands;
  - the tongue bands below the grid-point height;
  - the Step-1 spin-up check.

## Step 4 — Runs (background scripts, low effort; about 1.5 h)

- 3 variants (FA, FAS, FAG) × 3 parts (gentle / north / south) × 8 f_ref = **72 runs**, tag `fa1`.
- They start once the five-part runs (o5) finish, so we don't compete for the machines.
- They run on vierzack03/04/05/07 with the existing launcher pattern (`week3/chain_orient_parts_v2.sh` style).
- Parts are merged per hole (`merge_part_runs.py`), then scored. The new tag avoids stale `.done` skips.

## Step 5 — Results and records (main session)

- **Results:** in plain words, per band at depth, with figures.
- **Records:**
  - a new book chapter 13 ("The air from above"): the station test, the cause, the test and its result;
  - a decision-log entry (22 Oct: limitation reopened; the rule; the outcome);
  - memory updates (known limitations, a new free-air note).
- **If FAC is adopted:** Step 5 (the Bayesian calibration) runs on it.
- **The mass-balance stage, later and in its own plan:**
  - free air for the mass balance too, either as a free-air `.mdi` or by moving CentralEurope onto the per-node reader;
  - recalibration of all CentralEurope glaciers;
  - the GCM correction re-anchored (the GMIP4 `*_era5corrected` inputs must be re-corrected to the new reference);
  - whether glacier tongues need a cold-pool correction (Step 1.5);
  - the short precipitation investigation (above), because a new air changes the snow/rain split and so c_prec.

## Verification

- **Step 0:** each job's stored request matches the table; both files open and cover 1940-01 to 2025-12 with 18 levels
  and both variables.
- **Step 1:** the download matches Lander's file; the station rescore reproduces the station test to about 0.1–0.2 °C
  (expected, since the levels differ); the spin-up criterion is checked.
- **Step 2:** byte-identical off; mass balance identical on; the logged band air equals the table; one max-effort review.
- **Steps 3–4:** the rule commit hash is earlier than the first `fa1` log; all 48 runs finish; the merged runs have 1135
  points below 15 m, like o4.

## Glossary

| term | meaning |
|---|---|
| free air | the open air at a given height, read from ERA5's pressure levels away from the ground |
| 2 m air | ERA5's air 2 m above its own smooth model ground |
| pressure level / hPa | ERA5 stores the atmosphere on fixed pressures (hectopascal); 700 hPa is about 3000 m, 550 hPa about 5000 m |
| geopotential height | the height of a pressure level, from ERA5's geopotential divided by 9.80665 |
| grid point / node | the ERA5 0.25° cell a glacier takes its climate from; its height is `hclim` |
| dtdz | GloGEM's monthly lapse rate (12 values from ERA5 pressure levels, 1995–2022) |
| CDS | Copernicus Climate Data Store, ECMWF's download service; its API lets scripts fetch finished requests |
| FA / FAC | free air alone / free air with the surface cooling |
| Gilbert relation | surface = 0.95 × air − 3.4 °C, measured at Col du Dôme (Gilbert et al. 2014) |
| FRHS | the current central thermal variant (no T_off, ERA5 lapse rate, Herron–Langway meltwater, snow refreezing S) |
| o4 | the current scored run set: FRHS with three orientation parts per band |
| f_ref | share of meltwater that refreezes and heats; an "effective" parameter |
| error on left-out massifs | each massif predicted with the f_ref chosen on the others (M1) |
| zone | 3000–3800 m, the firn–no-firn transition |
| T_off | the mass-balance temperature correction per glacier; kept out of the thermal model |
| Herron–Langway (H-L), `hl_tk` | the firn densification scheme and the temperature it uses |
| gate / harness | checks that a switch changes nothing when off, and computes what it should when on |
| MB | mass balance |
| c_prec | the mass-balance calibration's precipitation factor per glacier |
| LAPrec | a station-based Alpine precipitation grid (MeteoSwiss), monthly |
| GLAMOS | Glacier Monitoring Switzerland, the source of the measured winter snow |
| GMIP4 | the GlacierMIP4 projections, whose climate inputs are already corrected to ERA5 |
