"""Which measured accumulations the perfect-burial test uses (fixed 21 Oct 2026 before its runs).

Input: the literature compilation of 21 Oct (week3/burial_research/measured_accumulation_by_entity.csv: one value per
scored entity, the mean over its holes with a value; sources in measured_accumulation_by_site.csv).
Used: entities with a measured long-term net accumulation above zero. Not used (the model's b_ref stays):
  - borrowed or basin values (Taconnaz, Jungfraufirn above 3415 m): not measured at the holes;
  - the Mont Blanc summit: two short records that conflict;
  - Hohsaas, Alphubel and Corvatsch: one recent year (2024/25, net ablation at most holes), not the long-term
    burial the reference balance stands for (a firn-extent question, covered by the observed-firn test);
  - entities without a measurement.
Writes accumulation.csv (entity, acc_central where used, acc_low, acc_high, confidence, use, flags) beside the input.

Usage: python scripts/icetemp_analysis/burial/select_accumulation.py [RESEARCH_DIR]
"""
import sys
from pathlib import Path

import pandas as pd

D = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/'
                                                         'burial_research')
OWN_ONE_YEAR = ('Alphubel South@', 'Hohsaasgletscher@', 'Vadret dal Corvatsch@')


def reason(r):
    if pd.isna(r.acc_entity_mean):
        return 'no measurement'
    if 'analogue' in str(r['flags']):
        return 'borrowed or basin value'
    if r.entity.startswith('Mont Blanc@'):
        return 'conflicting short records (Mont Blanc summit)'
    if r.entity.startswith(OWN_ONE_YEAR) or 'ablation' in str(r['flags']) or r.acc_entity_mean <= 0:
        return 'one recent year only (2024/25), not the long-term burial'
    return 'used'


def main():
    A = pd.read_csv(D / 'measured_accumulation_by_entity.csv')
    A['use'] = [reason(r) for _, r in A.iterrows()]
    out = A.assign(acc_central=A.acc_entity_mean.where(A.use == 'used'), confidence=A.weakest_confidence)
    out[['entity', 'acc_central', 'acc_low', 'acc_high', 'confidence', 'use', 'flags']].to_csv(
        D / 'accumulation.csv', index=False)
    print(out.use.value_counts().to_string())


if __name__ == '__main__':
    main()
