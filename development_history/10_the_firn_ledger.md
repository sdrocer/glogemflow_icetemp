# 10. The firn ledger: how much firn the meltwater can really reach

> **In short.** The thermal model lets meltwater warm firn down to a fixed depth, the Herron &
> Langway close-off depth of 75–90 m in the zone. Real firn is far thinner, and changes from year to
> year. The firn ledger keeps each band's annual firn layers from the mass balance and stops the water
> at their base. Scored on 14 Oct it was not adopted (chapter 9).
>
> This chapter asks why. First (check 5): the ledger's rules reproduce Janosch's firn model (Beer et
> al. 2026) when both get the same GLAMOS balances. But fed GloGEM's own balance, the ledger holds
> 10–30 m w.e. of firn at Hohsaas, Corvatsch and Alphubel where the observations leave almost none.
> The ledger is right; its input is not.
>
> *Status (2026-10-19): check 5 done. Next: the perfect-firn test (would the ledger help if the
> firn lines were right?), the rerun with the agreed amendments, and a firn column that follows the
> ledger (3b-ii). Plan: `plans/2026-10-19_firn_ledger_plan.md`.*

## What we believed

- **The ledger** (GloGEM branch `firn-ledger`, switch `firnice_perm_limit = 'firn'`), for each band and
  hydrological year:
  - a positive balance adds a layer;
  - a negative balance removes its amount from the newest layers (mass-conserving);
  - layers older than 20 years become ice;
  - the meltwater reaches at most the ledger plus this winter's snow.
- **It uses only GloGEM's own band balance,** so it runs everywhere GloGEM runs, including the
  projections.
- **On 14 Oct** it improved only one zone site (Tête Rousse) and cooled six young firn-line bands
  (chapter 9). We suspected the input: GloGEM's band balance is too positive on small high glaciers
  (Corvatsch, Hohsaas), so the ledger is too thick there.

```
              band balance, per year            rules                    firn per band and year
 firn model   GLAMOS distributed grids  ──────► firn_change.pro ────────► the reference
 check 5a     GLAMOS band balances      ──────► the ledger's rules ─────► same as the reference?   (the rules)
 check 5b     GloGEM's band balances    ──────► the ledger's rules ─────► how far off, and where   (the input)
```

## What we measured

### Check 5: the ledger against the firn model (2026-10-19)

**How** (`scripts/icetemp_analysis/firn/ledger_vs_firnmodel.py`; outputs in `week2/firnmodel/check5/`):
- **The firn model's numbers, taken apart one step at a time:**
  - (A) the firn model;
  - (B) the same GLAMOS cells with the ledger's rules, one rule at a time;
  - (C) the ledger's rules on the GLAMOS *band* balance;
  - (D) the ledger GloGEM keeps, fed its own balance (the FRHL run of 14 Oct).
- **The data:** the 9 glaciers with GLAMOS grids, on the RGI7 10 m bands with at least 3 GLAMOS cells,
  1972–2025.
- **The check that the comparison itself is right:** the rule code, set to the firn model's rules,
  reproduces the Python port of `firn_change.pro` exactly on all 9 glaciers. The firn differs by
  0.0 m w.e., and no cell gets a different firn/no-firn verdict.

![The ledger against the firn model in the borehole bands](figures/ch10_check5.svg)

**5a: the rules are right.** In the borehole bands, the ledger's rules on the GLAMOS balance follow the
firn model:

| firn in the borehole bands, 2000s (m w.e.) | firn model (A) | ledger rules, GLAMOS cells (B) | ledger rules, GLAMOS bands (C') | GloGEM's ledger (D) |
|---|---|---|---|---|
| Hohsaas 3275–3345 m | 0.2 | 0.5 | 0.3 | **24.5** |
| Corvatsch 3295–3345 m | 0.0 | 0.0 | 0.0 | **23.9** |
| Alphubel 3755–3825 m | 2.6 | 5.1 | 5.0 | **14.1** |
| Gr. Aletsch 3415–3635 m | 33.5 | 33.5 | 34.5 | 40.9 |

- **The rule differences are small:**
  - the mass-conserving removal (a decision of 14 Oct) keeps up to 2.5 m w.e. more (Alphubel,
    2000s);
  - the firn model's running-sum gate never changes anything;
  - working on band means instead of cells adds at most 1.0 m w.e. (Aletsch, 2000s).
- **The 20-year ageing is the same convention** in both (ages 0–20 are kept). The "21 layers" noted on
  14 Oct is not a difference.
- **The two-years-in-a-row rule** only decides when a cell counts as firn on a map; it does not
  change the firn's thickness, which is what the ledger uses. So nothing changes in GloGEM.
  - Counting firn from the first year, and the mass-conserving removal, put the firn line lower than
    the firn model's: by up to 175 m on Aletsch in the 2020s, mostly by tens of metres.
- **The start differs by design.** The firn model starts from no firn in 1970, the ledger from a
  steady state (here: of the first ten observed years). The two meet by the 1990s.

**5b: the input is wrong.** Fed GloGEM's band balance, the same rules give far too much firn on the
small high glaciers:
- **Hohsaas and Corvatsch:** 25–39 m w.e. through the 1980s and 1990s, and still 6–12 m w.e. in the
  2020s. The observations leave almost none after 2000.
- **Alphubel:** too thick until about 2015.
- **Gr. Aletsch:** much closer, but 12–18 m w.e. too thick before 2000 and 4 m w.e. too thin in the 2020s.
- **The firn line** from GloGEM's ledger lies lower than the GLAMOS-driven one, by 100–200 m on
  Hohsaas and Corvatsch in the 2000s–2010s. In the 2020s it lies higher on the large glaciers
  (Aletsch +157 m, Morteratsch +133 m).

**What it means:**
- The ledger's water limit is wrong at the small zone glaciers because GloGEM's balance builds firn
  that is not there, not because of the ledger's rules.
- The snowline recalibration did not fix that balance (chapter 9).
- The open question is whether the ledger would help if the firn were right. That is the
  perfect-firn test, next.

## What is next

1. **The perfect-firn test.** The same thermal model and ledger, but on the GLAMOS glaciers the
   ledger is fed the observed band balance:
   - FRHL-o: the thickness only;
   - FRHL-O: the thickness and the firn mask.

   If FRHL-O gains at least 2 robust zone sites over FRH (on Hohsaas, Corvatsch, Aletsch, Alphubel),
   investing in the snow and firn lines pays off.
2. **The rerun of the ledger** on the current code: the 360-day thermal year, the 3000–3800 m zone,
   and the refined cooling guard (mature firn only).
3. **3b-ii:** the firn density, and with it conduction and heat capacity, follows the ledger instead
   of one fixed profile.

## Sources

- Plan: `plans/2026-10-19_firn_ledger_plan.md`; the ledger's own plan `plans/2026-10-14_firn_ledger_plan.md`
- Scripts: `scripts/icetemp_analysis/firn/` (`firn_change_port.py`, `band_compare.py`,
  `ledger_vs_firnmodel.py`); `scripts/icetemp_analysis/firn_ledger.py` (the Python reference)
- GLAMOS distributed balances: `/scratch_net/vierzack06_fourth/mhuss/MassBalance_all/<glacier>/massbalance/results/`
- Literature: Beer et al. (2026, *The Cryosphere* 20, 5453–5473); Herron & Langway (1980)

## Glossary

| term | meaning |
|---|---|
| firn ledger | per elevation band, the annual firn layers built from the band balance; layers older than 20 years become ice; the meltwater reaches at most their base |
| firn model | Matthias Huss's `firn_change.pro` (Beer et al. 2026): the same bookkeeping per 25 m cell, driven by the GLAMOS distributed mass balance |
| GLAMOS | Glacier Monitoring Switzerland: measured mass balances, here as yearly grids |
| band balance | the year's mass balance of one 10 m elevation band [m w.e.] |
| m w.e. | metres of water equivalent |
| mass-conserving removal | a negative year removes exactly its amount, newest layers first (the firn model removes more when several thin layers lie on top; chapter 9) |
| running-sum gate | the firn model removes firn only where the running sum of the balance is positive |
| steady state | the starting ledger: one layer of a reference balance for each of the 21 years before the start |
| firn line | the lowest band where at least half the area has firn |
| borehole bands | the elevation bands of the boreholes on each glacier |
| zone | 3000–3800 m, the band of the firn–no-firn transition |
| FRH, FRHL | the central thermal variant (no T_off, ERA5 lapse rate, Herron & Langway meltwater); FRHL: with the ledger |
| FRHL-o, FRHL-O | the ledger fed the observed balance: o the thickness only, O also the firn mask |
| 3b-ii | the firn density follows the ledger instead of one fixed profile |
| robust zone site | a borehole site in the zone whose temperature class is right across the f_ref values |
