# Decision log

Every choice made on the borehole data that shaped the model or its calibration, with who made it, whether a rule
fixed before the runs existed, and whether the decision followed it. A referee should be able to see, for each
step, what was decided before and what after the data were looked at.

Kept since 13 Oct 2026 (the fair-test pre-registration); earlier changes are described in chapters 1–7.

| date | decision | by | rule fixed before? | followed the rule? | where |
|---|---|---|---|---|---|
| 13 Oct | Fair test of the thermal forcing: each forcing chooses its own f_ref leaving one massif out; scorer committed before the first run | Janosch, plan | yes (`plans/2026-10-13_fair_test_preregistration.md`, icetemp `fc469b2`) | — | ch. 8 |
| 14 Oct | FRH (no T_off in the thermal inputs, ERA5's own lapse rate) becomes the central forcing; H kept as a full sensitivity | rule A | yes | yes, but not robust (without Bernina H wins by 0.51 °C) | ch. 8 |
| 16 Oct | The thermal model stays on the monthly GloGEM (daily mass balance kept for later) | Janosch | no | — | ch. 9 |
| 19 Oct | The monthly snow-line recalibration is not used in the thermal runs | rule, Janosch agreed | yes | yes | ch. 9 |
| 19 Oct | Gain rules count sites, not bands; comparisons per band at depth, never site means | Janosch, before any result of the batch was read | amendment during the batch, before reading | — | ch. 10 |
| 20 Oct | The firn ledger is not adopted | stage-2 rule | yes | yes | ch. 10, 11 |
| 20 Oct | **S (meltwater refreezing in the snow on ice bands) adopted**; the central variant becomes FRHS | Janosch | yes (stage-2 rule) | **no**: the rule rejected S because one measured frozen bed (Tête Rousse) was lost; adopted because S is a real process and the remaining errors have other causes | ch. 11 |
| 20 Oct | W (station lapse rates) and B (full refreezing model) not adopted | stage rules | yes | yes | ch. 11 |
| 20 Oct | Comparisons pooled per borehole and calendar year (1135 points below 15 m) | plan | yes (plan of 20 Oct) | — | ch. 11 |
| 20 Oct | Test 1c (past mass-balance correction) reported, not adopted | plan | no rule; diagnostic | — | ch. 11 |
| 20 Oct | Step 5 plan revised after the review panel: error model chosen by evidence, coded site checks, nested benchmark against simple rules, f_ref lower limit 0.03, strain heating off, Chessjen sealed | Janosch | — | — | plan of 20 Oct |
| 21 Oct | Orientation adoption rule **B** (shade error inside sites at least halved, also without Taconnaz and Lys; no harm overall as for S) chosen over the investigation's stricter rule A; committed before the Stage-1 runs | Janosch | yes (icetemp `7bcbb8f`, `ab5ded4`) | — | ch. 12 |
| 21 Oct | Shared model bands take the point-weighted factor; three entities the band model cannot test for exposure leave the within-site measures | after a max-effort review, before any Stage-1 run was read | yes (`ab5ded4`) | — | ch. 12 |
| 21 Oct | f_ref per elevation class not pursued (worse on left-out massifs, 2.42 → 3.00 °C) | Janosch | no | — | plan of 20 Oct (revised 21 Oct) |
| 21 Oct | Stage 1: orientation **not met** by rule B (shade error −1.37 → −0.98, needed −0.68) | rule B | yes | yes | ch. 12 |
| 21 Oct | Perfect-burial test: the Mont Blanc / Monte Rosa split does **not** close (rule fixed before the runs, icetemp `990cd82`); over-burial reported as a known model error, no fix in v1.0 | rule; Janosch | yes | yes | ch. 12 |
| 21 Oct | Re-check of orientation on the measured burial: **not met** again (−1.47 → −1.04) | rule B | yes | yes | ch. 12 |
| 21 Oct | **Orientation kept in the model**, as three parts per band (gentle, north, south), three runs per glacier | Janosch | yes (rule B) | **no**: rule B failed twice; kept because it is physical and the map must show north against south faces. The paper claims the overall gain and only a weak within-site confirmation | ch. 12 |
| 21 Oct | The ERA5 air (about 2 °C too cold at all high stations, mostly in winter) and the missing snow-surface cooling stay as known limitations; no correction, no sensitivity run; the coupled fix is future work | Janosch | no | — | ch. 12 |
| 21 Oct | Tête Rousse: the drilled-thickness run dropped (the columns already reach the drilled depth); the colder spin-up tested instead, no effect, reported | Janosch | no rule; diagnostic | — | ch. 12 |
| 22 Oct | **The air reopened:** after the product search (ERA5's own free air +0.6 °C in winter at the high stations, today's air −5.2 °C), the free air is tested in the thermal model, with and without Gilbert's snow-surface cooling; mass balance unchanged. Replaces the 21 Oct "no correction, no sensitivity run" | Janosch | — | — | ch. 13 |
| 22 Oct | Free-air rule: FAG (free air + surface 0.95 × air − 1.4) is adopted if, against o4, the error on left-out massifs is not worse by more than 0.2 °C overall, above 3800 m and in the zone, regime classes lose at most one site net, and frozen beds (found minus false) drop by at most one; FA and FAS (−3.4, the measured skin) reported. −1.4 chosen because Gilbert's own firn model needed +2 K on the skin | Janosch | yes (icetemp `c66b30f`, before any `fa1` run) | — | ch. 13 |
| 22 Oct | Free-air rule, point 3 amended after the max-effort review, still before any `fa1` run: sites from the site table instead of the name before '@' (which made Grenzgletscher one group of 33 bands and Taconnaz two sites); the review also found the free air read one month late (fixed, GloGEM `d783f84`) | review; Claude | amendment before any run | — | ch. 13 |

## Not yet decided (Step 5)

- The error model (chosen by evidence among the listed forms, step 3.1 of the Step 5 plan).
- The site rule and site checks, the benchmark and its failure rule (step 3), all to be committed before any fit.
