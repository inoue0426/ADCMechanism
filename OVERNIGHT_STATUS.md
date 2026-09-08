# ADCMechanism Status

## Current Decision

**PIVOT-WEAK**

The matched-context pivot analysis finds descriptive within-context payload differences, but the current data cannot identify a defensible mechanism-level effect. The project should not proceed to a more complex mechanism-aware model on this dataset.

## What Was Tested

Commands:

```bash
python3 scripts/audit_adcdb_activity.py
python3 scripts/identifiability_study.py
python3 scripts/matched_context_pivot.py
```

Data:

- Local `data/raw/` was inspected; no replacement dataset was downloaded.
- ADCdb detail HTML pages: 435.
- Parsed ADCdb cell-line activity rows: 1,043.
- Usable deterministic molar IC50/LC50/EC50/GI50 rows at 10 nM: 705.
- Collapsed ADC/cell-line/antigen/payload observations: 624.
- Major-mechanism modeling subset excluding `Other`: 609 observations.
- Matched exact major-mechanism antigen x cell-line contexts: 24.
- Matched payload-context units after collapsing repeated payload observations: 58.
- Input observations represented by those matched units: 113.

## Strongest Result

Payload-held-out categorical baselines did not show a robust mechanism gain:

- `bio` = antigen + cell line: AP 0.8341, ROC-AUC 0.6805, MCC 0.2788.
- `bio_plus_mechanism`: AP 0.8294, ROC-AUC 0.6819, MCC 0.2159.
- `bio_plus_payload`: AP 0.8422, ROC-AUC 0.6931, MCC 0.2562.
- `bio_plus_payload_plus_mechanism`: AP 0.8446, ROC-AUC 0.6949, MCC 0.1852.

Payload leave-one-out agreed: `bio` AP 0.8333 vs `bio_plus_mechanism` AP 0.8203.

The final matched-context pivot found one descriptive pattern:

- DNA-damaging vs Topo-I inhibitor in 7 exact contexts: mean payload-unit response difference 0.3571, with 4/0/3 DNA-higher/Topo-higher/equal contexts.
- However, Topo-I is represented by only `PAY0ZVBAI` in this comparison.

## Strongest Confound Or Failure

- Payload identity deterministically maps to mechanism in this table: Cramer's V = 1.000.
- Mechanism is strongly associated with antigen, cell line, and cell-line/antigen context: Cramer's V = 0.828, 0.877, and 0.942.
- Only 24 crossed major-mechanism cell-line/antigen contexts exist, covering 113 observations.
- Only 58 payload-context units survive exact matched-context restriction, 13.2% of the major-mechanism payload-context units.
- Zero pairwise matched contexts have >=2 independent payload identities in both mechanism arms.
- Holding out `PAY0ZVBAI` eliminates the DNA-damaging vs Topo-I comparison entirely.
- True mechanism labels did not beat payload-level shuffled mechanism labels under the payload-group split: true AP 0.8294; shuffled mean AP 0.8261; shuffled p95 AP 0.8497.
- A context-fixed residual permutation was nominal for the all-mechanism residual range (observed 0.1874; p=0.0343), but pairwise mechanism permutation tests were not compelling: DNA vs microtubule p=1.0000, DNA vs Topo-I p=0.1266, Topo-I vs microtubule p=0.3294.

## Highest-Value Next Step

Stop #44 as a mechanism-aware modeling project on this dataset. If the idea is pursued further, the next experiment is data acquisition/curation: create exact antigen x cell-line contexts measured under at least two independent payloads per mechanism, across several contexts per mechanism pair, with harmonized endpoints and units.

## Main Outputs

- `results/identifiability_study.md`
- `results/cheap_model_results.csv`
- `results/shuffled_mechanism_control.csv`
- `results/matched_context_mechanism_effect_summary.csv`
- `results/matched_context_mechanism_effects.csv`
- `results/matched_context_pivot.md`
- `results/matched_context_pairwise_payload_unit_summary.csv`
- `results/matched_context_permutation_tests.csv`
- `results/matched_context_bootstrap_ci.csv`
- `results/matched_context_leave_one_payload_out_summary.csv`
- `data/processed/adc_response_10nM_context.csv`
