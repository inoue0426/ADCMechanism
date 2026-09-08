# ADC Mechanism Identifiability Study

## Scope

- Raw data were read from local `data/raw/`; no replacement dataset was downloaded.
- No neural model or complex architecture was trained.
- Parsed ADCdb cell-line activity rows available: 1,043.
- Usable deterministic molar potency rows at 10 nM: 705.
- Collapsed modeling observations: 624.
- Collapsed response distribution: positives=437, negatives=187, prevalence=0.700.
- Duplicate collapse conflicts at 10 nM: 11 context groups.

Usable rows require potency endpoint (IC50/LC50/EC50/GI50), molar unit convertible to nM, and a determinate 10 nM label. Mass units were excluded because molecular-weight conversion is not available in these local rows.

## Raw Data Inventory

| top_level | files | total_bytes | suffixes |
| --- | --- | --- | --- |
| adcdb_idrblab | 441 | 28649281 | .csv:1, .html:439, .txt:1 |
| adcnet_github | 4 | 219477 | .json:1, .xlsx:3 |
| annotationdb_adcdb | 2 | 67589 | .json:2 |
| digitaltumors_adcdb_pipeline | 3 | 6785 | .md:2, .yaml:1 |

## ADCNet Excel Files

| path | rows | columns | label_columns | label_counts |
| --- | --- | --- | --- | --- |
| data/raw/adcnet_github/data.xlsx | 435 | 15 | label（10nm）; label（100nm）; label（1nm）; label（1000nm） | label（10nm）:0=164,1=271; label（100nm）:0=154,1=281; label（1nm）:0=171,1=264; label（1000nm）:0=151,1=284 |
| data/raw/adcnet_github/files_data.xlsx | 435 | 15 | label（10nm）; label（100nm）; label（1nm）; label（1000nm） | label（10nm）:0=164,1=271; label（100nm）:0=154,1=281; label（1nm）:0=171,1=264; label（1000nm）:0=151,1=284 |
| data/raw/adcnet_github/t_data.xlsx | 50 | 15 | label（100nm）; label（1000nm） | label（100nm）:0=16,1=34; label（1000nm）:0=15,1=35 |

Local interpretation: `data.xlsx` and `files_data.xlsx` are byte-identical ADC-level ADCNet tables. They include threshold labels but no cell-line activity rows. `t_data.xlsx` is smaller and lacks a 10 nM label, so it is not the primary dataset for this question.

## Endpoint And Label Accounting

| reason | rows |
| --- | --- |
| ambiguous_censored_value | 13 |
| non_molar_or_missing_value | 246 |
| non_potency_endpoint | 79 |

## Per-Mechanism Audit

| payload_mechanism | observations | positives_10nM | prevalence_10nM | unique_payloads | effective_payload_n | unique_antigens | unique_cell_lines | top_payload | top_payload_share | top_antigen | top_antigen_share |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | 166 | 130 | 0.783 | 8 | 5.02 | 9 | 117 | PAY0IJSCK | 30.7% | Hepatocyte growth factor receptor (MET) | 30.7% |
| Other | 15 | 9 | 0.600 | 3 | 2.15 | 4 | 14 | PAY0VZNCV | 73.3% | Tumor necrosis factor receptor superfamily member 17 (TNFRSF17) | 73.3% |
| Topo-I inhibitor | 78 | 42 | 0.538 | 6 | 3.26 | 7 | 37 | PAY0ZVBAI | 62.8% | Receptor tyrosine-protein kinase erbB-2 (HER2) | 26.9% |
| microtubule-disrupting | 365 | 256 | 0.701 | 8 | 4.44 | 20 | 142 | PAY0FSXOW | 44.9% | Mast/stem cell growth factor receptor Kit (KIT) | 18.6% |

## Mechanism Associations

| association | right_unique_values | cramers_v |
| --- | --- | --- |
| mechanism_payload | 25 | 1.000 |
| mechanism_antigen | 28 | 0.828 |
| mechanism_cell_line | 256 | 0.877 |
| mechanism_disease | 93 | 0.694 |
| mechanism_cell_antigen_context | 382 | 0.942 |

Cramer's V near 1.0 indicates strong association. Payload identity almost deterministically implies mechanism, so seen-payload random splits cannot establish mechanism-level generalization.

## Cross-Mechanism Biological Overlap

### cell_line

| context | mechanism_count | unique_contexts |
| --- | --- | --- |
| cell_line | 1 | 209 |
| cell_line | 2 | 40 |
| cell_line | 3 | 7 |

### antigen_name

| context | mechanism_count | unique_contexts |
| --- | --- | --- |
| antigen_name | 1 | 20 |
| antigen_name | 2 | 5 |
| antigen_name | 3 | 2 |
| antigen_name | 4 | 1 |

### cell_antigen_context

| context | mechanism_count | unique_contexts |
| --- | --- | --- |
| cell_antigen_context | 1 | 357 |
| cell_antigen_context | 2 | 23 |
| cell_antigen_context | 3 | 2 |

Pairwise overlaps:

| context | mechanism_a | mechanism_b | shared_contexts |
| --- | --- | --- | --- |
| cell_line | DNA-damaging | microtubule-disrupting | 36 |
| cell_line | Topo-I inhibitor | microtubule-disrupting | 12 |
| cell_line | DNA-damaging | Topo-I inhibitor | 11 |
| cell_line | Other | microtubule-disrupting | 2 |
| cell_line | DNA-damaging | Other | 0 |
| cell_line | Other | Topo-I inhibitor | 0 |
| antigen_name | DNA-damaging | microtubule-disrupting | 4 |
| antigen_name | Topo-I inhibitor | microtubule-disrupting | 4 |
| antigen_name | DNA-damaging | Other | 3 |
| antigen_name | DNA-damaging | Topo-I inhibitor | 3 |
| antigen_name | Other | microtubule-disrupting | 2 |
| antigen_name | Other | Topo-I inhibitor | 1 |
| cell_antigen_context | DNA-damaging | microtubule-disrupting | 16 |
| cell_antigen_context | DNA-damaging | Topo-I inhibitor | 7 |
| cell_antigen_context | Topo-I inhibitor | microtubule-disrupting | 5 |
| cell_antigen_context | Other | microtubule-disrupting | 1 |
| cell_antigen_context | DNA-damaging | Other | 0 |
| cell_antigen_context | Other | Topo-I inhibitor | 0 |

Top crossed cell-line/antigen contexts:

| context | mechanism_count | observations | unique_payloads | has_label_variation | mechanisms |
| --- | --- | --- | --- | --- | --- |
| SK-OV-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 3 | 9 | 4 | 1 | DNA-damaging; Topo-I inhibitor; microtubule-disrupting |
| MDA-MB-231 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 3 | 7 | 4 | 0 | DNA-damaging; Topo-I inhibitor; microtubule-disrupting |
| SK-BR-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 11 | 4 | 0 | DNA-damaging; microtubule-disrupting |
| NCI-N87 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 10 | 4 | 1 | DNA-damaging; microtubule-disrupting |
| MCF-7 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 9 | 4 | 0 | Topo-I inhibitor; microtubule-disrupting |
| BT-474 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 6 | 4 | 1 | Other; microtubule-disrupting |
| A375.S2 cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| CVCL_0479 || Cadherin-6 (CDH6) | 2 | 5 | 2 | 0 | Topo-I inhibitor; microtubule-disrupting |
| Calu-6 cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| Hs 700T cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| MDA-MB-468 cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| NCI-H1703 cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| PA-1 cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| Raji cells || CD276 antigen (CD276) | 2 | 5 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| MDA-MB-231 cells || Tyrosine-protein kinase receptor UFO (AXL) | 2 | 3 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| Raji cells || B-cell receptor CD22 (CD22) | 2 | 3 | 2 | 1 | DNA-damaging; Topo-I inhibitor |
| Ramos cells || B-cell receptor CD22 (CD22) | 2 | 3 | 2 | 1 | DNA-damaging; Topo-I inhibitor |
| Reh cells || B-cell receptor CD22 (CD22) | 2 | 3 | 2 | 1 | DNA-damaging; Topo-I inhibitor |
| U-87MG cells || Epidermal growth factor receptor (EGFR) | 2 | 3 | 2 | 0 | DNA-damaging; microtubule-disrupting |
| CVCL_0033 || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 2 | 2 | 0 | Topo-I inhibitor; microtubule-disrupting |

## Matched-Context Mechanism Effects

This uses only observations where the same cell-line/antigen context has more than one major mechanism. It is a low-power stress test for mechanism effects after holding both measured cell line and ADC antigen fixed.

| mechanism_a | mechanism_b | crossed_contexts | contexts_with_label_variation | observations_in_pair_contexts | mean_prevalence_diff_a_minus_b | weighted_prevalence_diff_a_minus_b | a_higher_contexts | b_higher_contexts | equal_contexts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | microtubule-disrupting | 16 | 2 | 76 | 0.0069 | -0.0906 | 1 | 1 | 14 |
| DNA-damaging | Topo-I inhibitor | 7 | 4 | 23 | 0.3571 | 0.4130 | 4 | 0 | 3 |
| Topo-I inhibitor | microtubule-disrupting | 5 | 1 | 30 | -0.2000 | -0.2667 | 0 | 1 | 4 |

Largest matched-context differences:

| cell_antigen_context | mechanism_a | mechanism_b | n_a | n_b | prevalence_a | prevalence_b | prevalence_diff_a_minus_b | payloads_a | payloads_b |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SK-OV-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | Topo-I inhibitor | 1 | 4 | 1.0000 | 0.0000 | 1.0000 | PAY0QKJME | PAY0ZVBAI |
| Raji cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 2 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Ramos cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 2 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Reh cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 2 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Daudi cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0RZSDY | PAY0ZVBAI |
| MDA-MB-231 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | Topo-I inhibitor | 1 | 4 | 0.0000 | 0.0000 | 0.0000 | PAY0QKJME | PAY0ZVBAI |
| RS4 11 cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0RZSDY | PAY0ZVBAI |
| SW620 cells || Epidermal growth factor receptor (EGFR) | DNA-damaging | microtubule-disrupting | 1 | 1 | 1.0000 | 0.0000 | 1.0000 | PAY0AKDAM | PAY0FSXOW |
| NCI-N87 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | microtubule-disrupting | 1 | 9 | 0.0000 | 0.8889 | -0.8889 | PAY0QKJME | PAY0FSXOW;PAY0JCIBW;PAY0VZVPQ |
| A375.S2 cells || CD276 antigen (CD276) | DNA-damaging | microtubule-disrupting | 1 | 4 | 1.0000 | 1.0000 | 0.0000 | PAY0SOSMQ | PAY0FSXOW |
| Calu-6 cells || CD276 antigen (CD276) | DNA-damaging | microtubule-disrupting | 1 | 4 | 1.0000 | 1.0000 | 0.0000 | PAY0SOSMQ | PAY0FSXOW |
| HCT 15 cells || Epidermal growth factor receptor (EGFR) | DNA-damaging | microtubule-disrupting | 1 | 1 | 0.0000 | 0.0000 | 0.0000 | PAY0AKDAM | PAY0FSXOW |
| Hs 700T cells || CD276 antigen (CD276) | DNA-damaging | microtubule-disrupting | 1 | 4 | 1.0000 | 1.0000 | 0.0000 | PAY0SOSMQ | PAY0FSXOW |
| LoVo cells || Epidermal growth factor receptor (EGFR) | DNA-damaging | microtubule-disrupting | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0AKDAM | PAY0FSXOW |
| MDA-MB-231 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | microtubule-disrupting | 1 | 2 | 0.0000 | 0.0000 | 0.0000 | PAY0QKJME | PAY0FSXOW;PAY0JCIBW |

## Phase 2 Decision

- Decision before modeling: **GO_FOR_CHEAP_FALSIFICATION**

- major mechanisms each have at least 5 payloads (minimum 6).
- major mechanisms each have at least 50 observations (minimum 78).
- there are 24 crossed cell-line/antigen contexts covering 113 observations.
- 6 crossed contexts contain both response classes.
- payload identity almost perfectly determines mechanism; random-split mechanism gains would not be identifiable.

## Cheap Modeling Falsification

Model: categorical Naive Bayes with Laplace smoothing. Splits are deterministic. The central split is grouped by payload identity, so test payloads are unseen during training.

| split | model | average_precision | roc_auc | balanced_accuracy | mcc | log_loss | brier |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cell_antigen_context_group | antigen_only | 0.8180 | 0.7187 | 0.6533 | 0.2941 | 0.5603 | 0.1841 |
| cell_antigen_context_group | bio | 0.7921 | 0.6719 | 0.6537 | 0.2853 | 0.6244 | 0.2108 |
| cell_antigen_context_group | bio_plus_mechanism | 0.7929 | 0.6754 | 0.6576 | 0.2926 | 0.6308 | 0.2123 |
| cell_antigen_context_group | bio_plus_payload | 0.8204 | 0.7073 | 0.6547 | 0.2997 | 0.6459 | 0.2109 |
| cell_antigen_context_group | bio_plus_payload_plus_mechanism | 0.8179 | 0.7032 | 0.6456 | 0.2741 | 0.6710 | 0.2170 |
| cell_antigen_context_group | cell_only | 0.6665 | 0.4748 | 0.4395 | -0.1218 | 0.6836 | 0.2391 |
| cell_antigen_context_group | global | 0.6762 | 0.4506 | 0.5000 |  | 0.6108 | 0.2098 |
| cell_antigen_context_group | mechanism_only | 0.7526 | 0.5566 | 0.5504 | 0.1378 | 0.6033 | 0.2069 |
| cell_antigen_context_group | payload_only | 0.8244 | 0.6618 | 0.5870 | 0.1841 | 0.5709 | 0.1945 |
| cell_line_group | antigen_only | 0.8110 | 0.6911 | 0.6569 | 0.3132 | 0.5836 | 0.1943 |
| cell_line_group | bio | 0.8124 | 0.6928 | 0.6282 | 0.2345 | 0.6047 | 0.2039 |
| cell_line_group | bio_plus_mechanism | 0.8001 | 0.6782 | 0.6254 | 0.2323 | 0.6134 | 0.2064 |
| cell_line_group | bio_plus_payload | 0.8277 | 0.6964 | 0.6536 | 0.2869 | 0.6307 | 0.2120 |
| cell_line_group | bio_plus_payload_plus_mechanism | 0.8228 | 0.6925 | 0.6377 | 0.2614 | 0.6554 | 0.2179 |
| cell_line_group | cell_only | 0.6491 | 0.4503 | 0.5000 |  | 0.6330 | 0.2203 |
| cell_line_group | global | 0.6491 | 0.4503 | 0.5000 |  | 0.6107 | 0.2098 |
| cell_line_group | mechanism_only | 0.7362 | 0.5493 | 0.5428 | 0.0941 | 0.6122 | 0.2109 |
| cell_line_group | payload_only | 0.8192 | 0.6549 | 0.5582 | 0.1167 | 0.5800 | 0.1987 |
| payload_group | antigen_only | 0.7892 | 0.6186 | 0.6187 | 0.2175 | 0.6231 | 0.2114 |
| payload_group | bio | 0.8341 | 0.6805 | 0.6467 | 0.2788 | 0.6112 | 0.2077 |
| payload_group | bio_plus_mechanism | 0.8294 | 0.6819 | 0.6107 | 0.2159 | 0.6014 | 0.2024 |
| payload_group | bio_plus_payload | 0.8422 | 0.6931 | 0.6401 | 0.2562 | 0.7054 | 0.2487 |
| payload_group | bio_plus_payload_plus_mechanism | 0.8446 | 0.6949 | 0.6009 | 0.1852 | 0.6700 | 0.2337 |
| payload_group | cell_only | 0.7745 | 0.6207 | 0.6203 | 0.2202 | 0.5900 | 0.2000 |
| payload_group | global | 0.6192 | 0.3353 | 0.5000 |  | 0.6276 | 0.2168 |
| payload_group | mechanism_only | 0.6815 | 0.3977 | 0.4879 | -0.0253 | 0.6421 | 0.2224 |
| payload_group | payload_only | 0.6835 | 0.3840 | 0.5000 |  | 0.6816 | 0.2442 |
| payload_leave_one_out | antigen_only | 0.7885 | 0.6160 | 0.6203 | 0.2227 | 0.6226 | 0.2114 |
| payload_leave_one_out | bio | 0.8333 | 0.6788 | 0.6186 | 0.2287 | 0.6116 | 0.2080 |
| payload_leave_one_out | bio_plus_mechanism | 0.8203 | 0.6666 | 0.6133 | 0.2239 | 0.6170 | 0.2083 |
| payload_leave_one_out | bio_plus_payload | 0.8361 | 0.6847 | 0.6357 | 0.2488 | 0.7103 | 0.2501 |
| payload_leave_one_out | bio_plus_payload_plus_mechanism | 0.8259 | 0.6741 | 0.6062 | 0.1951 | 0.6935 | 0.2431 |
| payload_leave_one_out | cell_only | 0.7872 | 0.6359 | 0.6261 | 0.2307 | 0.5869 | 0.1987 |
| payload_leave_one_out | global | 0.6204 | 0.3163 | 0.5000 |  | 0.6215 | 0.2143 |
| payload_leave_one_out | mechanism_only | 0.6393 | 0.3694 | 0.4038 | -0.2027 | 0.6505 | 0.2257 |
| payload_leave_one_out | payload_only | 0.6709 | 0.3380 | 0.5000 |  | 0.6805 | 0.2437 |
| random_row | antigen_only | 0.8435 | 0.7452 | 0.6691 | 0.3361 | 0.5267 | 0.1739 |
| random_row | bio | 0.9037 | 0.8051 | 0.7159 | 0.4248 | 0.4718 | 0.1559 |
| random_row | bio_plus_mechanism | 0.9013 | 0.8016 | 0.7090 | 0.4252 | 0.4770 | 0.1575 |
| random_row | bio_plus_payload | 0.9021 | 0.8051 | 0.7366 | 0.4733 | 0.4966 | 0.1632 |
| random_row | bio_plus_payload_plus_mechanism | 0.8986 | 0.8005 | 0.7266 | 0.4451 | 0.5169 | 0.1677 |
| random_row | cell_only | 0.8605 | 0.7507 | 0.6570 | 0.2870 | 0.5187 | 0.1720 |
| random_row | global | 0.6266 | 0.4973 | 0.5000 |  | 0.6085 | 0.2089 |
| random_row | mechanism_only | 0.6996 | 0.5631 | 0.5504 | 0.1378 | 0.6003 | 0.2052 |
| random_row | payload_only | 0.8154 | 0.6791 | 0.6135 | 0.2498 | 0.5609 | 0.1891 |

Shuffle control: mechanism labels were permuted at payload level and evaluated under the same payload-group split.

- Shuffle repeats: 50
- Shuffled AP mean: 0.8261
- Shuffled AP min/max: 0.7779 / 0.8709

## Current Decision

- **PIVOT**

- payload-group AP gain for bio+mechanism over bio: -0.0047.
- payload-group AP gain for bio+payload+mechanism over bio+payload: 0.0024.
- payload-group MCC gain for bio+mechanism over bio: -0.0629.
- true bio+mechanism AP 0.8294; shuffled-mechanism mean 0.8261, p95 0.8497.
- true mechanism does not beat the shuffled-mechanism control.
- matched-context stress test has 7/28 pairwise contexts with label variation.

## Pivot Recommendation

The defensible reformulation is narrower: test whether coarse payload mechanism can act as a low-dimensional prior for unseen payloads in ADCdb-like metadata, not whether it independently explains ADC response beyond payload identity and biological context.

## Highest-Value Next Step

Stop model escalation. The next useful step is curation, not architecture: validate cell-line naming/synonyms, endpoint aggregation, and mass-unit conversion assumptions to see whether more genuinely matched cell-line/antigen contexts can be recovered from the existing local ADCdb pages.
