# ADC Mechanism-Aware Response Prediction: Data Audit Pilot

## Scope

- Model training: not performed.
- ADCdb detail HTML pages parsed: 435
- Cell-line activity rows parsed: 1,043
- Cell-line IC50/EC50/GI50 potency rows used for audit tables: 964
- Unique ADCs with parsed potency activity: 173
- Unique payloads in potency rows: 51
- Unique antigens in potency rows: 36
- Unique cell lines in potency rows: 319
- Processed activity CSV: `data/processed/adcdb_cell_line_activity.csv`

Source scope: ADC IDs were taken from the public ADCNet `data.xlsx` file, then each ID was resolved against the public ADCdb detail page. This audits the ADCNet/ADCdb subset rather than a complete ADCdb dump.

Candidate data sources checked:

- Official ADCdb detail pages: chosen because they expose ADC, payload, antigen, payload target, and cell-line potency rows.
- ADCNet GitHub data.xlsx: used only as the ADCdb ID seed; the table has 10 nM labels but does not retain cell-line activity rows.
- AnnotationDB ADC API / digitaltumors ADCdb pipeline: useful for ADC/component metadata, but the tested individual JSON did not include cell-line activity rows.
- Zenodo ADCdb mirror: identified, but direct record/file downloads returned HTTP 504 during this run.

## 1. Sample Counts Per Payload Mechanism

| payload_mechanism | activity_rows | nM_rows | molar_unit_rows | unique_adcs | unique_payloads | unique_antigens | unique_cell_lines | top_payload_share | top_payload | top_units |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | 256 | 122 | 185 | 27 | 13 | 10 | 154 | 71/256 (27.7%) | PAY0AKDAM | nM:122, ng/mL:62, pM:60, ug/mL:8 |
| Other | 81 | 14 | 16 | 38 | 23 | 6 | 18 | 11/81 (13.6%) | PAY0VZNCV | ug/mL:63, nM:14, ng/mL:2, uM:2 |
| Topo-I inhibitor | 149 | 118 | 119 | 19 | 6 | 7 | 43 | 80/149 (53.7%) | PAY0ZVBAI | nM:118, ng/mL:27, ug/mL:3, uM:1 |
| microtubule-disrupting | 478 | 352 | 398 | 89 | 9 | 27 | 178 | 227/478 (47.5%) | PAY0FSXOW | nM:352, ng/mL:37, ug/mL:29, pM:27 |

## 2. Unique Payloads Per Mechanism

| payload_mechanism | unique_payloads | activity_rows | top_payload | top_payload_rows | top_payload_fraction |
| --- | --- | --- | --- | --- | --- |
| DNA-damaging | 13 | 256 | PAY0AKDAM | 71 | 27.7% |
| Other | 23 | 81 | PAY0VZNCV | 11 | 13.6% |
| Topo-I inhibitor | 6 | 149 | PAY0ZVBAI | 80 | 53.7% |
| microtubule-disrupting | 9 | 478 | PAY0FSXOW | 227 | 47.5% |

Payload-mechanism consistency diagnostic:

| unique_payloads | payloads_seen_in_one_mechanism | fraction_payloads_seen_in_one_mechanism | payloads_seen_in_multiple_mechanisms |
| --- | --- | --- | --- |
| 51 | 51 | 1.000 | 0 |

Interpretation: payloads are expected to map mostly to one mechanism. The stronger leakage risk is when a mechanism class is represented by only one payload or is dominated by one payload.

Payload target mapping used for mechanism assignment:

| payload_target | assigned_mechanism | activity_rows |
| --- | --- | --- |
| Microtubule (MT) | microtubule-disrupting | 478 |
| Human deoxyribonucleic acid (hDNA) | DNA-damaging | 159 |
| DNA topoisomerase 1 (TOP1) | Topo-I inhibitor | 149 |
| Human Deoxyribonucleic acid (hDNA) | DNA-damaging | 82 |
| Glucocorticoid receptor (NR3C1) | Other | 64 |
| DNA topoisomerase 2-alpha (TOP2A) | DNA-damaging | 15 |
| RNA polymerase II | Other | 11 |
| Stimulator of interferon genes protein (STING1) | Other | 3 |
| Phospholipid hydroperoxide glutathione peroxidase (GPX4) | Other | 2 |
| Ribosome (RB) | Other | 1 |

## 3. Antigen And Cell Line Counts Per Mechanism

| payload_mechanism | unique_antigens | unique_cell_lines | unique_adcs |
| --- | --- | --- | --- |
| DNA-damaging | 10 | 154 | 27 |
| Other | 6 | 18 | 38 |
| Topo-I inhibitor | 7 | 43 | 19 |
| microtubule-disrupting | 27 | 178 | 89 |

## 4. Class Imbalance

- Max/min sample-count ratio across non-empty mechanism classes: 5.90x

## 5. Cross-Mechanism Cell Line Coverage

- Cell lines measured in >1 mechanism class: 62 / 319
- Activity rows belonging to those cross-mechanism cell lines: 414 / 964

| mechanism_classes_per_cell_line | unique_cell_lines |
| --- | --- |
| 1 | 257 |
| 2 | 51 |
| 3 | 10 |
| 4 | 1 |

Pairwise shared cell-line counts:

| mechanism_a | mechanism_b | shared_cell_lines |
| --- | --- | --- |
| DNA-damaging | microtubule-disrupting | 47 |
| Topo-I inhibitor | microtubule-disrupting | 19 |
| DNA-damaging | Topo-I inhibitor | 14 |
| Other | microtubule-disrupting | 4 |
| Other | Topo-I inhibitor | 2 |
| DNA-damaging | Other | 1 |

Top cross-mechanism cell-line examples:

| cell_line | mechanism_count | activity_rows | mechanisms |
| --- | --- | --- | --- |
| MDA-MB-468 cells | 4 | 16 | DNA-damaging, Other, Topo-I inhibitor, microtubule-disrupting |
| SK-BR-3 cells | 3 | 33 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| NCI-N87 cells | 3 | 28 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| MDA-MB-231 cells | 3 | 21 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| Raji cells | 3 | 20 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| Ramos cells | 3 | 20 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| BT-474 cells | 3 | 19 | Other, Topo-I inhibitor, microtubule-disrupting |
| MCF-7 cells | 3 | 18 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| SK-OV-3 cells | 3 | 18 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| Daudi cells | 3 | 14 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| SK-MES-1 cells | 3 | 4 | DNA-damaging, Topo-I inhibitor, microtubule-disrupting |
| T-47D cells | 2 | 14 | DNA-damaging, microtubule-disrupting |
| Reh cells | 2 | 10 | DNA-damaging, Topo-I inhibitor |
| Karpas-299 cells | 2 | 9 | DNA-damaging, microtubule-disrupting |
| HL-60 cells | 2 | 8 | DNA-damaging, microtubule-disrupting |

## Decision

- 判定: **検証可能**
- 主因: 主要な即時停止理由なし

Next comparison design:

- Use identical train/validation/test splits for baseline and mechanism-conditioned models.
- Split at least by cell line, and additionally report held-out-payload sensitivity if sample size allows.
- Baseline: antigen/cell-line/delivery features without mechanism input.
- Mechanism-conditioned: same features plus payload mechanism indicator or mechanism-specific head.
- Primary comparison: PR-AUC/MCC/balanced accuracy on the 10 nM response label, stratified by mechanism and by held-out cell lines.
