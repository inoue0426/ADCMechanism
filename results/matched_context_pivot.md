# Matched-Context Pivot Analysis

## Scope

- Input: existing local `data/processed/adc_response_10nM_context.csv`.
- Raw data were not redownloaded or replaced.
- No predictive model, neural head, foundation model, or hyperparameter search was run.
- Central analysis excludes `Other` and uses the three interpretable major mechanisms: DNA-damaging, Topo-I inhibitor, and microtubule-disrupting.
- Unit of inference: one payload identity within one exact antigen x cell-line context. Repeated ADC observations using the same payload in the same context are collapsed before comparison.

## Data Reduction

| scope | input_observations | positive_input_observations | negative_input_observations | observation_prevalence | payload_context_units | mean_payload_unit_response_rate | unique_payloads | unique_mechanisms | unique_cell_antigen_contexts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_processed_observations | 624 | 437 | 187 | 0.7003 | 453 | 0.7165 | 25 | 4 | 382 |
| major_mechanisms_only | 609 | 428 | 181 | 0.7028 | 439 | 0.7199 | 22 | 3 | 369 |
| matched_major_exact_contexts | 113 | 80 | 33 | 0.7080 | 58 | 0.7131 | 13 | 3 | 24 |
| matched_major_survival_fraction_vs_major | 18.6% |  |  |  | 13.2% |  |  |  |  |

## Exact Matched-Context Coverage

- Exact crossed major-mechanism antigen x cell-line contexts: 24.
- Payload-context units in those contexts: 58.
- Input observations represented by those units: 113.
- Crossed contexts with any unit-level response variation: 6.
- Pairwise context comparisons where both mechanism arms contain >=2 independent payload identities: 0.
- Strong-overlap contexts, defined as >=3 mechanisms or >=4 payload-context units: 5.

| cell_antigen_context | mechanism_count | payload_context_units | input_observations | has_unit_response_variation | payload_units_by_mechanism | unit_response_by_mechanism | payloads_by_mechanism |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MDA-MB-231 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 3 | 4 | 7 | 0 | DNA-damaging=1; Topo-I inhibitor=1; microtubule-disrupting=2 | DNA-damaging=0.0000; Topo-I inhibitor=0.0000; microtubule-disrupting=0.0000 | DNA-damaging=PAY0QKJME; Topo-I inhibitor=PAY0ZVBAI; microtubule-disrupting=PAY0FSXOW,PAY0JCIBW |
| SK-OV-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 3 | 4 | 9 | 1 | DNA-damaging=1; Topo-I inhibitor=1; microtubule-disrupting=2 | DNA-damaging=1.0000; Topo-I inhibitor=0.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0QKJME; Topo-I inhibitor=PAY0ZVBAI; microtubule-disrupting=PAY0FSXOW,PAY0VZVPQ |
| MCF-7 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 4 | 9 | 0 | Topo-I inhibitor=1; microtubule-disrupting=3 | Topo-I inhibitor=0.0000; microtubule-disrupting=0.0000 | Topo-I inhibitor=PAY0ZVBAI; microtubule-disrupting=PAY0FSXOW,PAY0JCIBW,PAY0VZVPQ |
| NCI-N87 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 4 | 10 | 1 | DNA-damaging=1; microtubule-disrupting=3 | DNA-damaging=0.0000; microtubule-disrupting=0.9524 | DNA-damaging=PAY0QKJME; microtubule-disrupting=PAY0FSXOW,PAY0JCIBW,PAY0VZVPQ |
| SK-BR-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 4 | 11 | 0 | DNA-damaging=1; microtubule-disrupting=3 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0LZRXU; microtubule-disrupting=PAY0FSXOW,PAY0JCIBW,PAY0VZVPQ |
| A375.S2 cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| CVCL_0033 || Receptor tyrosine-protein kinase erbB-2 (HER2) | 2 | 2 | 2 | 0 | Topo-I inhibitor=1; microtubule-disrupting=1 | Topo-I inhibitor=1.0000; microtubule-disrupting=1.0000 | Topo-I inhibitor=PAY0WSZPJ; microtubule-disrupting=PAY0FSXOW |
| CVCL_0479 || Cadherin-6 (CDH6) | 2 | 2 | 5 | 0 | Topo-I inhibitor=1; microtubule-disrupting=1 | Topo-I inhibitor=1.0000; microtubule-disrupting=1.0000 | Topo-I inhibitor=PAY0UXDSB; microtubule-disrupting=PAY0GTSVM |
| Calu-6 cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| Daudi cells || B-cell receptor CD22 (CD22) | 2 | 2 | 2 | 0 | DNA-damaging=1; Topo-I inhibitor=1 | DNA-damaging=1.0000; Topo-I inhibitor=1.0000 | DNA-damaging=PAY0RZSDY; Topo-I inhibitor=PAY0ZVBAI |
| HCT 15 cells || Epidermal growth factor receptor (EGFR) | 2 | 2 | 2 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=0.0000; microtubule-disrupting=0.0000 | DNA-damaging=PAY0AKDAM; microtubule-disrupting=PAY0FSXOW |
| Hs 700T cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| LoVo cells || Epidermal growth factor receptor (EGFR) | 2 | 2 | 2 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0AKDAM; microtubule-disrupting=PAY0FSXOW |
| MDA-MB-231 cells || Tyrosine-protein kinase receptor UFO (AXL) | 2 | 2 | 3 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0XBVCC; microtubule-disrupting=PAY0FSXOW |
| MDA-MB-468 cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| NCI-H1703 cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| PA-1 cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| RS4 11 cells || B-cell receptor CD22 (CD22) | 2 | 2 | 2 | 0 | DNA-damaging=1; Topo-I inhibitor=1 | DNA-damaging=1.0000; Topo-I inhibitor=1.0000 | DNA-damaging=PAY0RZSDY; Topo-I inhibitor=PAY0ZVBAI |
| Raji cells || B-cell receptor CD22 (CD22) | 2 | 2 | 3 | 1 | DNA-damaging=1; Topo-I inhibitor=1 | DNA-damaging=1.0000; Topo-I inhibitor=0.5000 | DNA-damaging=PAY0RZSDY; Topo-I inhibitor=PAY0ZVBAI |
| Raji cells || CD276 antigen (CD276) | 2 | 2 | 5 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=0.0000; microtubule-disrupting=0.0000 | DNA-damaging=PAY0SOSMQ; microtubule-disrupting=PAY0FSXOW |
| Ramos cells || B-cell receptor CD22 (CD22) | 2 | 2 | 3 | 1 | DNA-damaging=1; Topo-I inhibitor=1 | DNA-damaging=1.0000; Topo-I inhibitor=0.5000 | DNA-damaging=PAY0RZSDY; Topo-I inhibitor=PAY0ZVBAI |
| Reh cells || B-cell receptor CD22 (CD22) | 2 | 2 | 3 | 1 | DNA-damaging=1; Topo-I inhibitor=1 | DNA-damaging=1.0000; Topo-I inhibitor=0.5000 | DNA-damaging=PAY0RZSDY; Topo-I inhibitor=PAY0ZVBAI |
| SW620 cells || Epidermal growth factor receptor (EGFR) | 2 | 2 | 2 | 1 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=0.0000 | DNA-damaging=PAY0AKDAM; microtubule-disrupting=PAY0FSXOW |
| U-87MG cells || Epidermal growth factor receptor (EGFR) | 2 | 2 | 3 | 0 | DNA-damaging=1; microtubule-disrupting=1 | DNA-damaging=1.0000; microtubule-disrupting=1.0000 | DNA-damaging=PAY0AKDAM; microtubule-disrupting=PAY0FSXOW |

## Within-Context Mechanism Effects

The response difference is computed from collapsed payload-context units, then averaged over exact contexts. Observation-weighted columns are reported as a sensitivity check, not as independent mechanism evidence.

| mechanism_a | mechanism_b | crossed_contexts | payload_context_units_in_pair_contexts | contexts_with_unit_response_difference | a_higher_contexts | b_higher_contexts | equal_contexts | mean_context_unit_response_diff_a_minus_b | observation_weighted_diff_a_minus_b | contexts_with_both_arms_2plus_payloads | unique_payloads_a_across_contexts | unique_payloads_b_across_contexts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | microtubule-disrupting | 16 | 38 | 2 | 1 | 1 | 14 | 0.0030 | -0.0906 | 0 | 5 | 3 |
| DNA-damaging | Topo-I inhibitor | 7 | 14 | 4 | 4 | 0 | 3 | 0.3571 | 0.4130 | 0 | 2 | 1 |
| Topo-I inhibitor | microtubule-disrupting | 5 | 14 | 1 | 0 | 1 | 4 | -0.2000 | -0.2667 | 0 | 3 | 4 |

Largest individual exact-context differences:

| cell_antigen_context | mechanism_a | mechanism_b | payload_units_a | payload_units_b | unit_response_a | unit_response_b | unit_response_diff_a_minus_b | payloads_a | payloads_b |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SK-OV-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 0.0000 | 1.0000 | PAY0QKJME | PAY0ZVBAI |
| SW620 cells || Epidermal growth factor receptor (EGFR) | DNA-damaging | microtubule-disrupting | 1 | 1 | 1.0000 | 0.0000 | 1.0000 | PAY0AKDAM | PAY0FSXOW |
| SK-OV-3 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | Topo-I inhibitor | microtubule-disrupting | 1 | 2 | 0.0000 | 1.0000 | -1.0000 | PAY0ZVBAI | PAY0FSXOW;PAY0VZVPQ |
| NCI-N87 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | microtubule-disrupting | 1 | 3 | 0.0000 | 0.9524 | -0.9524 | PAY0QKJME | PAY0FSXOW;PAY0JCIBW;PAY0VZVPQ |
| Raji cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Ramos cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Reh cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 0.5000 | 0.5000 | PAY0RZSDY | PAY0ZVBAI |
| Daudi cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0RZSDY | PAY0ZVBAI |
| MDA-MB-231 cells || Receptor tyrosine-protein kinase erbB-2 (HER2) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 0.0000 | 0.0000 | 0.0000 | PAY0QKJME | PAY0ZVBAI |
| RS4 11 cells || B-cell receptor CD22 (CD22) | DNA-damaging | Topo-I inhibitor | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0RZSDY | PAY0ZVBAI |
| A375.S2 cells || CD276 antigen (CD276) | DNA-damaging | microtubule-disrupting | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0SOSMQ | PAY0FSXOW |
| Calu-6 cells || CD276 antigen (CD276) | DNA-damaging | microtubule-disrupting | 1 | 1 | 1.0000 | 1.0000 | 0.0000 | PAY0SOSMQ | PAY0FSXOW |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

## Context-Fixed Residual Check

Residuals subtract each exact context's mean payload-unit response before summarizing by mechanism.

| payload_mechanism | payload_context_units | unique_payloads | unique_contexts | mean_context_fixed_residual | median_context_fixed_residual | min_context_fixed_residual | max_context_fixed_residual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | 21 | 6 | 21 | 0.0374 | 0.0000 | -0.7143 | 0.5000 |
| Topo-I inhibitor | 10 | 3 | 10 | -0.1500 | 0.0000 | -0.7500 | 0.0000 |
| microtubule-disrupting | 27 | 4 | 19 | 0.0265 | 0.0000 | -0.5000 | 0.2857 |

## Negative Controls

Permutation shuffles mechanism labels only within exact contexts, preserving context and mechanism arm sizes. The payload bootstrap resamples payload identities as clusters; degenerate rows indicate that one side has fewer than two payloads, so payload-level uncertainty is not estimable.

| test | mechanism_a | mechanism_b | observed_statistic | permutation_mean | permutation_p025 | permutation_p975 | p_two_sided_or_upper_tail |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all_mechanisms_context_fixed_residual_range |  |  | 0.1874 | 0.0940 | 0.0121 | 0.1874 | 0.0343 |
| pair_context_mean_unit_response_diff | DNA-damaging | microtubule-disrupting | 0.0030 | -0.0004 | -0.1220 | 0.0863 | 1.0000 |
| pair_context_mean_unit_response_diff | DNA-damaging | Topo-I inhibitor | 0.3571 | -0.0011 | -0.3571 | 0.3571 | 0.1266 |
| pair_context_mean_unit_response_diff | Topo-I inhibitor | microtubule-disrupting | -0.2000 | 0.0012 | -0.2000 | 0.1000 | 0.3294 |

| bootstrap | mechanism_a | mechanism_b | observed_statistic | bootstrap_p025 | bootstrap_p975 | fraction_gt_zero | payload_bootstrap_degenerate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| context_and_within_arm_payload_unit_bootstrap | DNA-damaging | microtubule-disrupting | 0.0030 | -0.1786 | 0.1875 | 0.4852 | 0 |
| payload_cluster_bootstrap | DNA-damaging | microtubule-disrupting | 0.0030 | -0.3333 | 0.1667 | 0.4907 | 0 |
| context_and_within_arm_payload_unit_bootstrap | DNA-damaging | Topo-I inhibitor | 0.3571 | 0.1429 | 0.6429 | 0.9968 | 0 |
| payload_cluster_bootstrap | DNA-damaging | Topo-I inhibitor | 0.3571 | 0.3000 | 0.5000 | 1.0000 | 1 |
| context_and_within_arm_payload_unit_bootstrap | Topo-I inhibitor | microtubule-disrupting | -0.2000 | -0.6000 | 0.0000 | 0.0000 | 0 |
| payload_cluster_bootstrap | Topo-I inhibitor | microtubule-disrupting | -0.2000 | -0.3333 | 0.0000 | 0.0000 | 0 |

Leave-one-payload-out stability:

| mechanism_a | mechanism_b | payloads_tested | unavailable_holdouts | payloads_making_statistic_unavailable | sign_change_holdouts | min_remaining_crossed_contexts | max_remaining_crossed_contexts | min_leave_one_payload_out_statistic | max_leave_one_payload_out_statistic |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | Topo-I inhibitor | 3 | 1 | PAY0ZVBAI | 0 | 0 | 5 | 0.3000 | 0.5000 |
| DNA-damaging | microtubule-disrupting | 8 | 0 |  | 2 | 4 | 16 | -0.2500 | 0.0769 |
| Topo-I inhibitor | microtubule-disrupting | 7 | 0 |  | 1 | 2 | 5 | -0.2500 | 0.0000 |

Strong-overlap subset summary:

| mechanism_a | mechanism_b | crossed_contexts | contexts_with_unit_response_difference | mean_context_unit_response_diff_a_minus_b | contexts_with_both_arms_2plus_payloads | unique_payloads_a_across_contexts | unique_payloads_b_across_contexts |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DNA-damaging | microtubule-disrupting | 4 | 1 | -0.2381 | 0 | 2 | 3 |
| Topo-I inhibitor | microtubule-disrupting | 3 | 1 | -0.3333 | 0 | 1 | 3 |
| DNA-damaging | Topo-I inhibitor | 2 | 1 | 0.5000 | 0 | 1 | 1 |

## Decision

**PIVOT-WEAK**

- DNA-damaging vs microtubule-disrupting: 16 exact contexts, mean unit-response diff 0.0030, 1/1/14 higher/lower/equal, both-arm >=2-payload contexts=0, permutation p=1.0000.
- DNA-damaging vs Topo-I inhibitor: 7 exact contexts, mean unit-response diff 0.3571, 4/0/3 higher/lower/equal, both-arm >=2-payload contexts=0, permutation p=0.1266.
- Topo-I inhibitor vs microtubule-disrupting: 5 exact contexts, mean unit-response diff -0.2000, 0/1/4 higher/lower/equal, both-arm >=2-payload contexts=0, permutation p=0.3294.
- No mechanism pair has replicated independent payload support within the same exact context (total both-arm >=2-payload pair-contexts=0).
- There is a descriptive within-context pattern, but it is not identifiable as a mechanism-level effect because one or both arms are driven by single payloads.

Interpretation: the strongest exact matched-context pattern is descriptive, not identifiable as a mechanism-level effect. The data can show payload-specific differences within some contexts, but they do not replicate across multiple independent payloads per mechanism inside the same contexts.

Additional crossed data required to resolve the pivot: the same antigen x cell-line contexts measured under at least two independent payloads per mechanism, across several contexts per mechanism pair, with harmonized endpoints and units. Without that, mechanism remains inseparable from payload identity in this dataset.
