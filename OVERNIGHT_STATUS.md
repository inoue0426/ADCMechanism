# ADCMechanism Status

## Current Decision

**PIVOT**

The local ADCdb/ADCNet data can support a cheap falsification test, but the current evidence does not justify building a larger mechanism-aware architecture.

## What Was Tested

Commands:

```bash
python3 scripts/audit_adcdb_activity.py
python3 scripts/identifiability_study.py
```

Data:

- Local `data/raw/` was inspected; no replacement dataset was downloaded.
- ADCdb detail HTML pages: 435.
- Parsed ADCdb cell-line activity rows: 1,043.
- Usable deterministic molar IC50/LC50/EC50/GI50 rows at 10 nM: 705.
- Collapsed ADC/cell-line/antigen/payload observations: 624.
- Major-mechanism modeling subset excluding `Other`: 609 observations.

## Strongest Result

Payload-held-out categorical baselines did not show a robust mechanism gain:

- `bio` = antigen + cell line: AP 0.8341, ROC-AUC 0.6805, MCC 0.2788.
- `bio_plus_mechanism`: AP 0.8294, ROC-AUC 0.6819, MCC 0.2159.
- `bio_plus_payload`: AP 0.8422, ROC-AUC 0.6931, MCC 0.2562.
- `bio_plus_payload_plus_mechanism`: AP 0.8446, ROC-AUC 0.6949, MCC 0.1852.

Payload leave-one-out agreed: `bio` AP 0.8333 vs `bio_plus_mechanism` AP 0.8203.

## Strongest Confound Or Failure

- Payload identity deterministically maps to mechanism in this table: Cramer's V = 1.000.
- Mechanism is strongly associated with antigen, cell line, and cell-line/antigen context: Cramer's V = 0.828, 0.877, and 0.942.
- Only 24 crossed major-mechanism cell-line/antigen contexts exist, covering 113 observations.
- Only 7/28 matched pairwise mechanism contexts have label variation.
- True mechanism labels did not beat payload-level shuffled mechanism labels under the payload-group split: true AP 0.8294; shuffled mean AP 0.8261; shuffled p95 AP 0.8497.

## Highest-Value Next Step

Do curation before any model escalation: validate cell-line naming/synonyms, endpoint aggregation, and mass-unit conversion assumptions to see whether additional genuinely matched cell-line/antigen contexts can be recovered from the existing local ADCdb pages.

## Main Outputs

- `results/identifiability_study.md`
- `results/cheap_model_results.csv`
- `results/shuffled_mechanism_control.csv`
- `results/matched_context_mechanism_effect_summary.csv`
- `results/matched_context_mechanism_effects.csv`
- `data/processed/adc_response_10nM_context.csv`
